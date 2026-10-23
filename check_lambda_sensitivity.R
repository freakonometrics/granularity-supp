
local({
  ROOT <- getwd()
  OUT <- file.path(ROOT,"checks")
  dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
  if (!requireNamespace("nloptr",quietly=TRUE))
    stop('Install nloptr: install.packages("nloptr")')
  SNAPSHOT <- file.path(ROOT,"engine_v3_snapshot.R")
  if (!file.exists(SNAPSHOT)) stop("Place engine_v3_snapshot.R next to this script")
  # Local R environment; no writes into the already running FULL session.
  eng <- new.env(parent=globalenv())
  sys.source(SNAPSHOT,envir=eng)
  stopifnot(exists("run_market_simulation_revised",envir=eng))
  engine_md5 <- unname(tools::md5sum(SNAPSHOT))
  T <- as.integer(Sys.getenv("CLIMATE_LAMBDA_T","120"))
  raw_l <- strsplit(Sys.getenv("CLIMATE_LAMBDA_VALUES","0.20,0.35,0.70"),",")[[1]]
  lambdas <- unique(as.numeric(trimws(raw_l)))
  stopifnot(length(lambdas)>=2L,all(is.finite(lambdas)),all(lambdas>0),all(lambdas<=1),
            is.finite(T),T>=20L)
  msg <- function(...) message(sprintf(...))
  beta <- 1; I <- 20L; R <- 3L
  SEED <- 123L; N_POP <- 500L
  set.seed(SEED)
  d <- matrix(rgamma(I*R,shape=rep(1,each=I)),I,R)
  omega <- d/rowSums(d)
  cost <- outer(as.vector(omega%*%c(.1,.6,.9)),10*c(.1,.6,.9)*c(.2,.5,.8))
  mass <- N_POP*omega
  p_pool <- sum(cost*mass)/sum(mass)+1.5/beta
  if(abs(p_pool-3.22980467)>1e-5) stop("Seed/cost generator mismatch with manuscript")
  msg("R2.3 pooled Nash price = %.10f; guard horizon = %d",p_pool,T)
  firms <- c("A","B","C")

  base_args <- list(
    seed=SEED,I=I,N_POP=N_POP,ALPHA=10,ETA_CAT=1,
    RISK_PROB=c(.1,.6,.9),RISK_COST_PROB=c(.2,.5,.8),
    BETA_PRICE=beta,V_INSURER=c(A=0,B=0,C=0),
    ARCHETYPES=c(A="national",B="communal",C="communal_risk"),
    INITIAL_QUOTES=c(A=p_pool,B=p_pool,C=p_pool),
    EPS_SHARE=.01,EPS_CUMULATIVE=.25,DAMPING_NU=.7,
    PRICE_LB=.01,PRICE_UB=Inf,N_STARTS=5L,
    quiet=TRUE
  )
  # Cache keyed by engine checksum, lambda, scenario and horizon.
  run_one <- function(scenario,lam,tmax,guard=FALSE) {
    id <- sprintf("%s_lambda_%s_T%d_%s",scenario,
      gsub("\\.","p",sprintf("%.2f",lam)),tmax,substr(engine_md5,1,10))
    f <- file.path(OUT,paste0(id,".rds"))
    if(file.exists(f)) {msg("CACHE %s",id);return(readRDS(f))}
    msg("RUN %s",id)
    args <- c(base_args,list(scenario=scenario,LAMBDA_REALLOC=lam,
      T_MAX=as.integer(tmax),TOL_CONVERGENCE=if(guard)0 else 1e-8,
      CUMULATIVE_PRESENCE=guard,
      FEASIBILITY_GUARD=guard,GUARD_TOL=1e-8,
      JOINT_REPAIR=guard,REPAIR_NSTARTS=4L,REPAIR_MAXEVAL=600L))
    result <- do.call(eng$run_market_simulation_revised,args)
    if(guard)eng$check_market_invariants(result,tol=1e-8)
    saveRDS(result,f)
    result
  }
  # Convert one engine output to insurer-level trajectories, using BOOK pbar.
  summarize_by_period <- function(z,scenario,lam) {
    ans <- list()
    per <- sort(unique(z$period))
    for(p in per) {
      x <- z[z$period==p,,drop=FALSE]
      hidx <- x$risk_class==max(x$risk_class)
      for(f in firms) {
        E <- x[[paste0(f,"_exposure")]]
        q <- x[[paste0(f,"_prix")]]
        pb <- x[[paste0(f,"_pbar")]]
        claims <- sum(E*x$pure_cost)
        ans[[length(ans)+1L]] <- data.frame(
          scenario=scenario,lambda=lam,iteration=p,insurer=f,
          quote_mean=sum(E*q)/sum(E),
          book_mean=sum(E*pb)/sum(E),
          highrisk_exposure=sum(E[hidx]),
          highrisk_share=sum(E[hidx])/sum(x$segment_mass[hidx]),
          book_profit=sum(E*(pb-x$pure_cost)),
          highrisk_book_profit=sum(E[hidx]*(pb[hidx]-x$pure_cost[hidx])))
      }
    }
    do.call(rbind,ans)
  }
  # Independent stock-flow replay. The UNCONSTRAINED quote path is lambda-
  # invariant analytically: the inherited part is additive in the current
  # optimization. After its last computed quote, freeze at the converged quote
  # and continue to the common horizon T. This avoids rerunning SQP 3 times.
  replay_unconstrained <- function(base,lam,T) {
    stop_at <- max(base$period)
    S <- matrix(1/3,nrow=I*R,ncol=3,dimnames=list(NULL,firms))
    PB <- matrix(p_pool,nrow=I*R,ncol=3,dimnames=list(NULL,firms))
    # Order all rows consistently by commune then risk class.
    extract <- function(t) {
      z <- base[base$period==min(t,stop_at),,drop=FALSE]
      z <- z[order(z$commune_id,z$risk_class),,drop=FALSE]
      if(nrow(z)!=I*R)stop("Unexpected number of cells in unconstrained path")
      z
    }
    records <- vector("list",T+1L)
    for(t in 0:T) {
      rowlist <- vector("list",length(firms))
      z <- extract(t)
      q <- cbind(A=z$A_prix,B=z$B_prix,C=z$C_prix)
      if(t>0L) {
        util <- -beta*q
        util <- util-apply(util,1,max)
        target <- exp(util)
        target <- target/rowSums(target)
        SN <- (1-lam)*S+lam*target
        PBN <- ((1-lam)*S*PB+lam*target*q)/SN
        S <- SN;PB <- PBN
      }
      for(j in seq_along(firms)) {
        f <- firms[j];E <- z$segment_mass*S[,j]
        ind <- z$risk_class==max(z$risk_class)
        rowlist[[j]] <- data.frame(scenario="unconstrained_replay",
          lambda=lam,iteration=t,insurer=f,
          quote_mean=sum(E*q[,j])/sum(E),
          book_mean=sum(E*PB[,j])/sum(E),
          highrisk_exposure=sum(E[ind]),
          highrisk_share=sum(E[ind])/sum(z$segment_mass[ind]),
          book_profit=sum(E*(PB[,j]-z$pure_cost)),
          highrisk_book_profit=sum(E[ind]*(PB[ind,j]-z$pure_cost[ind])))
      }
      records[[t+1L]] <- do.call(rbind,rowlist)
    }
    do.call(rbind,records)
  }
  # For the reference unconstrained path use a short engine run, which
  # ordinarily reaches its stopping criterion near iteration 36. This runs
  # independently of the simultaneously executing 200-period figure job.
  base <- run_one("unconstrained",.35,80L,guard=FALSE)
  msg("Unconstrained quote stopping iteration: %d",max(base$period))
  all_records <- list()
  for(lam in lambdas)
    all_records[[length(all_records)+1L]] <- replay_unconstrained(base,lam,T)
  write.csv(do.call(rbind,all_records),file.path(OUT,"R23_unconstrained_lambda_paths.csv"),
            row.names=FALSE)
  # Truly constrained runs must be computed separately because lambda enters
  # the individual feasible set and the coordinated acceptance/projection step.
  guarded_diagnostics <- list()
  for(lam in lambdas) {
    msg("Starting safeguarded COMMUNE path: lambda %.2f, T=%d",lam,T)
    z <- tryCatch(run_one("commune",lam,T,guard=TRUE),error=function(e) e)
    if(inherits(z,"error")) {
      warning(sprintf("Failed lambda %.2f: %s",lam,conditionMessage(z)))
      writeLines(conditionMessage(z),file.path(OUT,sprintf("ERROR_lambda_%0.2f.txt",lam)))
      next
    }
    all_records[[length(all_records)+1L]] <- summarize_by_period(z,"commune",lam)
    cv <- eng$convergence_diagnostics(z)
    gd <- eng$guard_diagnostics(z)
    sl <- eng$constraint_diagnostics(z)
    last <- cv[nrow(cv),,drop=FALSE]
    guarded_diagnostics[[length(guarded_diagnostics)+1L]] <- data.frame(
      lambda=lam,last_iteration=max(z$period),
      quote_residual=last$max_price_change,
      share_residual=last$max_share_change,
      total_backtracks=sum(gd$step_scale<1,na.rm=TRUE),
      total_joint_projections=sum(gd$joint_repair,na.rm=TRUE),
      min_accepted_slack=if(nrow(sl))min(sl$slack)else NA_real_,
      n_violations=if(nrow(sl))sum(sl$slack< -1e-8)else 0L,
      stationary=(is.finite(last$max_price_change) && is.finite(last$max_share_change) &&
                    last$max_price_change<1e-7 && last$max_share_change<1e-7))
    utils::write.csv(gd,file.path(OUT,sprintf("R23_repairs_lambda_%0.2f.csv",lam)),row.names=FALSE)
  }
  all_paths <- do.call(rbind,all_records)
  write.csv(all_paths,file.path(OUT,"R23_all_paths.csv"),row.names=FALSE)
  terminal <- all_paths[!duplicated(all_paths[c("scenario","lambda","insurer")],
                                  fromLast=TRUE),,drop=FALSE]
  write.csv(terminal,file.path(OUT,"R23_lambda_terminal_prices_shares_profits.csv"),
            row.names=FALSE)
  if(length(guarded_diagnostics))
    write.csv(do.call(rbind,guarded_diagnostics),
              file.path(OUT,"R23_guarded_commune_diagnostics.csv"),row.names=FALSE)
  writeLines(c(
    "R2.3 lambda sensitivity: engine snapshot kept separate from ongoing FULL job",
    paste("Engine MD5:",engine_md5),paste("R:",R.version.string),
    paste("nloptr version:",as.character(utils::packageVersion("nloptr"))),
    sprintf("Seed=%d I=%d beta=1 nu=0.7 epsilon=0.01 epsilon_cum=0.25",SEED,I),
    paste("Guard horizon:",T,"lambdas:",paste(lambdas,collapse=",")),
    "Unconstrained: EXACT stock-flow replay of the lambda-independent quote path",
    "after baseline quote stopping, quoted prices are held at their final value.",
    "Constrained commune: separate solver runs; fixed horizon, guard+projection.",
    "Do not describe guarded terminal states as Nash or stationary unless certified."
  ),file.path(OUT,"R23_run_notes.txt"))
  message("R2.3 DONE. Results in ",normalizePath(OUT))
  print(terminal[terminal$insurer=="C",],row.names=FALSE)
})

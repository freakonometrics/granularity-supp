# =====================================================================
# R2.1 — independent scalar best-response validation, no SQP / no nloptr.
# Reproduces exactly the 20-commune, seed-123 cost and Dirichlet generator.
# Architecture: A national, B communal, C commune x risk cell.
# Output is SEPARATE from the ongoing 200-iteration FULL calculations.
# Run: source("check_R21_best_responses.R")
# =====================================================================
local({
  options(stringsAsFactors = FALSE)
  ROOT <- getwd()
  OUT <- file.path(ROOT, "checks_R21")
  dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
  sink_path <- file.path(OUT, "run_notes.txt")
  SEED <- 123L; I <- 20L; R <- 3L; N_POP <- 500
  beta <- 1; nu <- .7
  N_GRID_ITER <- as.integer(Sys.getenv("CLIMATE_R21_GRID_ITER", "241"))
  N_GRID_FINAL <- as.integer(Sys.getenv("CLIMATE_R21_GRID_FINAL", "1501"))
  MAX_ITER <- as.integer(Sys.getenv("CLIMATE_R21_MAX_ITER", "160"))
  stopifnot(N_GRID_ITER >= 80L, N_GRID_FINAL >= 400L, MAX_ITER > 1L)
  set.seed(SEED)
  raw <- matrix(rgamma(I*R, shape=rep(1,each=I)), I, R)
  omega <- raw / rowSums(raw)
  u <- as.vector(omega %*% c(.1,.6,.9))
  mass <- N_POP * omega
  cost <- outer(u, 10*c(.1,.6,.9)*c(.2,.5,.8))
  bar_cost <- sum(mass*cost) / sum(mass)
  p_pool <- bar_cost + 1.5 / beta
  firms <- c("A","B","C")
  if(abs(p_pool-3.22980467)>1e-5) stop("Seed/cost generator mismatch with manuscript")
  message(sprintf("R2.1: pooled quote %.10f; 20 communes; 60 risk cells",p_pool))

  # Fully stable cellwise logit for any two rival prices.
  sigma_at <- function(q, riv1, riv2) {
    # With v_i=0, sigma_i(q)=1/(1+[exp(-beta*p_j)+exp(-beta*p_k)] exp(beta*q)).
    # p_j,p_k finite and positive; clipping prevents under/overflow in search.
    k <- exp(-beta * riv1) + exp(-beta * riv2)
    plogis(-beta*q - log(k))
  }
  target_shares <- function(P) {
    denom <- exp(-beta*P$A)+exp(-beta*P$B)+exp(-beta*P$C)
    setNames(lapply(firms,function(nm) exp(-beta*P[[nm]])/denom),firms)
  }
  # Scalar optimization for a fixed tariff group. Group cells are weighted
  # with fixed population mass, not an endogenously changing insurer mix.
  scalar_br <- function(nm,cell_indices,P,n_grid) {
    other <- setdiff(firms,nm)
    rc <- cost[cell_indices]; w <- mass[cell_indices]
    r1 <- P[[other[1]]][cell_indices]; r2 <- P[[other[2]]][cell_indices]
    rivals <- c(r1,r2)
    # Finite upper search bound with a NEGATIVE derivative at and beyond it:
    # q>=max(cost)+4/beta and q>=max(rival)+4/beta imply small sigma and
    # a negative derivative for each cell (hence no distant optimum).
    lo <- .01
    hi <- max(10, max(rc)+4/beta, max(rivals)+4/beta)
    probability <- function(q) sigma_at(q,r1,r2)
    f <- function(q) {
      sh <- probability(q)
      sum(w*sh*(q-rc))
    }
    derivative <- function(q) {
      sh <- probability(q)
      sum(w*sh*(1-beta*(1-sh)*(q-rc)))
    }
    gx <- seq(lo,hi,length.out=n_grid)
    gd <- vapply(gx,derivative,numeric(1))
    crossings <- which(gd[-length(gd)] * gd[-1] <= 0)
    roots <- if(length(crossings)) vapply(crossings,function(j)
      uniroot(derivative,lower=gx[j],upper=gx[j+1],tol=1e-11)$root,
      numeric(1)) else numeric()
    candidates <- unique(c(lo,hi,roots))
    values <- vapply(candidates,f,numeric(1))
    j <- which.max(values)
    list(q=candidates[j],profit=values[j],n_roots=length(roots),
         deriv_upper=derivative(hi),search_upper=hi,
         min_grid_deriv=min(gd),max_grid_deriv=max(gd),
         endpoint_winner=j<=2L)
  }
  global_br <- function(P,n_grid) {
    mats <- list(A=matrix(NA_real_,I,R), B=matrix(NA_real_,I,R),
                 C=matrix(NA_real_,I,R))
    details <- list()
    put <- function(nm,key,ix,z) {
      mats[[nm]][ix] <<- z$q
      details[[length(details)+1L]] <<- data.frame(
        insurer=nm,tariff_cell=key,BR_quote=z$q,BR_group_profit=z$profit,
        n_stationary_roots=z$n_roots,upper_search_price=z$search_upper,
        derivative_at_upper=z$deriv_upper,endpoint_winner=z$endpoint_winner)
    }
    ix <- seq_len(I*R)
    put("A","national",ix,scalar_br("A",ix,P,n_grid))
    for(m in seq_len(I)) {
      ix <- m + (seq_len(R)-1L)*I # row m, 3 columns
      put("B",sprintf("commune_%02d",m),ix,scalar_br("B",ix,P,n_grid))
    }
    for(m in seq_len(I))for(r in seq_len(R)) {
      ix <- m+(r-1L)*I
      put("C",sprintf("m%02d_r%d",m,r),ix,scalar_br("C",ix,P,n_grid))
    }
    list(prices=mats,details=do.call(rbind,details))
  }
  payoffs <- function(P) {
    shares <- target_shares(P)
    vapply(firms,function(nm)sum(mass*shares[[nm]]*(P[[nm]]-cost)),numeric(1))
  }
  highrisk <- function(P) {
    shares <- target_shares(P)
    vapply(firms,function(nm)sum(mass[,3]*shares[[nm]][,3])/sum(mass[,3]),numeric(1))
  }
  # The independently computed iterates need not coincide step-by-step with
  # the paper's SQP; the relevant test is their terminal price and regret.
  P <- setNames(lapply(firms,function(nm) matrix(p_pool,I,R)),firms)
  path <- list()
  for(t in seq_len(MAX_ITER)) {
    best <- global_br(P,N_GRID_ITER)
    Pnew <- setNames(lapply(firms,function(nm)
      (1-nu)*P[[nm]]+nu*best$prices[[nm]]),firms)
    residual <- max(vapply(firms,function(nm)
      max(abs(Pnew[[nm]]-P[[nm]])),numeric(1)))
    P <- Pnew
    h <- highrisk(P)
    path[[t]] <- data.frame(iteration=t,price_change=residual,
                           A_H=h['A'],B_H=h['B'],C_H=h['C'])
    if(t==1L || t%%5L==0L || residual < 1e-9)
      message(sprintf("R2.1 iter %d: price residual %.3g; C high-risk %.5f",
                      t,residual,h['C']))
    if(residual < 1e-9) break
  }
  # Use a FINER independent scan for all 81 tariff coefficients at the end.
  final <- global_br(P,N_GRID_FINAL)
  u_current <- payoffs(P)
  # Instead change one insurer at a time, holding others at the observed profile.
  u_br <- vapply(firms,function(nm) {
    unilateral <- P;unilateral[[nm]] <- final$prices[[nm]]
    payoffs(unilateral)[[nm]]
  },numeric(1))
  regret <- pmax(0,u_br-u_current)
  diff <- vapply(firms,function(nm)
    max(abs(final$prices[[nm]]-P[[nm]])),numeric(1))
  ff <- data.frame(insurer=firms,profit_at_profile=as.numeric(u_current),
                   profit_at_global_scalar_BR=as.numeric(u_br),
                   unilateral_gain=as.numeric(regret),
                   relative_gain=as.numeric(regret)/pmax(1,abs(u_current)),
                   max_price_BR_residual=as.numeric(diff))
  check <- data.frame(mean_cost=bar_cost, pooled_price=p_pool,
                      C_H_share=unname(highrisk(P)['C']),
                      C_H_exposure=unname(highrisk(P)['C'])*sum(mass[,3]),
                      final_max_price_change=path[[length(path)]][["price_change"]],
                      all_upper_derivatives_negative=all(final$details$derivative_at_upper<0),
                      no_search_endpoint_winners=!any(final$details$endpoint_winner),
                      MAX_ITER=MAX_ITER,N_GRID_FINAL=N_GRID_FINAL)
  terminal <- do.call(rbind,lapply(firms,function(nm) {
    d <- expand.grid(commune=seq_len(I),risk_class=seq_len(R))
    data.frame(insurer=nm,commune=d$commune,risk_class=d$risk_class,
      weight=as.vector(mass),cost=as.vector(cost),
      quote=as.vector(P[[nm]]),scalar_best_response=as.vector(final$prices[[nm]]),
      price_deviation=as.vector(final$prices[[nm]]-P[[nm]]))
  }))
  # Optional: test the ACTUAL SQP terminal quotes, when its RDS is present.
  # This check is stronger evidence for the paper than comparing only our
  # independently computed terminal allocation. It is SAFE while another
  # R process is writing caches: a failed/incomplete RDS is simply skipped.
  sqp_file <- Sys.getenv("CLIMATE_R21_SOURCE_RDS", "")
  if(!nzchar(sqp_file)) {
    d0 <- file.path(ROOT,"climate_extension_outputs","full","seed_123","runs")
    candidates <- if(dir.exists(d0))
      list.files(d0,pattern="^(account_unconstrained|hetero_from_pool)_.*[.]rds$",
                 full.names=TRUE) else character()
    if(length(candidates))
      sqp_file <- candidates[order(file.info(candidates)$mtime,decreasing=TRUE)][1]
  }
  if(nzchar(sqp_file) && file.exists(sqp_file)) {
    observed <- tryCatch(readRDS(sqp_file),error=function(e)NULL)
    if(!is.null(observed) &&
       all(c("period","commune_id","risk_class","segment_mass","pure_cost",
             "A_prix","B_prix","C_prix") %in% names(observed))) {
      z <- observed[observed$period==max(observed$period),,drop=FALSE]
      z <- z[order(z$commune_id,z$risk_class),,drop=FALSE]
      if(nrow(z)==I*R &&
         max(abs(z$segment_mass-as.vector(t(mass))))<1e-6 &&
         max(abs(z$pure_cost-as.vector(t(cost))))<1e-6) {
        Pobs <- lapply(firms,function(nm) {
          matrix(z[[paste0(nm,"_prix")]],nrow=I,ncol=R,byrow=TRUE)
        })
        names(Pobs) <- firms
        test <- global_br(Pobs,N_GRID_FINAL)
        obsp <- payoffs(Pobs)
        obsbrp <- vapply(firms,function(nm) {
          trial<-Pobs;trial[[nm]]<-test$prices[[nm]]
          payoffs(trial)[[nm]]
        },numeric(1))
        sqp <- data.frame(source_rds=basename(sqp_file),
          source_terminal_iteration=max(z$period),insurer=firms,
          profit_at_SQP=as.numeric(obsp),
          unilateral_best_profit=as.numeric(obsbrp),
          numerical_regret=pmax(0,as.numeric(obsbrp-obsp)),
          BR_quote_residual=vapply(firms,function(nm)
            max(abs(Pobs[[nm]]-test$prices[[nm]])),numeric(1)))
        write.csv(sqp,file.path(OUT,"R21_SQP_profile_regrets.csv"),row.names=FALSE)
        message("Also validated existing SQP RDS: ",basename(sqp_file))
      } else warning("SQP RDS uses a different market; ignored: ",sqp_file)
    } else warning("SQP RDS unreadable or schema mismatch: ",sqp_file)
  } else {
    message("No SQP RDS provided; independent static profile was validated.")
  }
  write.csv(do.call(rbind,path),file.path(OUT,"R21_convergence_path.csv"),row.names=FALSE)
  write.csv(ff,file.path(OUT,"R21_regrets_by_insurer.csv"),row.names=FALSE)
  write.csv(final$details,file.path(OUT,"R21_scalar_global_checks.csv"),row.names=FALSE)
  write.csv(terminal,file.path(OUT,"R21_terminal_prices.csv"),row.names=FALSE)
  write.csv(check,file.path(OUT,"R21_market_validation.csv"),row.names=FALSE)
  saveRDS(list(prices=P,regrets=ff,market=check),file.path(OUT,"R21_equilibrium_check.rds"))
  writeLines(c(paste("R version:",R.version.string),
    sprintf("Seed %d / communes %d / beta %.6f",SEED,I,beta),
    "Method: all scalar derivative sign changes on a finite grid + uniroot.",
    "The search upper bound lies in a region where all cell derivatives are negative.",
    "A finite grid does not formally exclude closely spaced missed roots.",
    "This is an independent numerical best-response check, not a formal proof of uniqueness.",
    paste("All roots scanned:",sum(final$details$n_stationary_roots)),
    paste("All endpoint winners:",sum(final$details$endpoint_winner))),sink_path)
  message("R2.1 DONE. Files in ",normalizePath(OUT))
  print(ff,row.names=FALSE)
  print(check,row.names=FALSE)
})

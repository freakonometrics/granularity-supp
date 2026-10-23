#| label: embedded-simulation-engine
#| echo: true
# ============================================================================
# Complete, visible simulation engine for the revised climate insurance paper.
# Standalone reimplementation of the stock-flow logit pricing model documented
# in supp-main.qmd and run_climate_extensions.R. No external source() or CSV.
# Requires: nloptr (SQP); optionally ggplot2 for visualization.
# Note: independently verify regeneration against the archived full outputs.
# ============================================================================

run_market_simulation_revised <- function(
  scenario=c("unconstrained","national","commune","risk_class"),
  seed=123L, ARCHETYPES=c(A="national",B="communal",C="communal_risk"),
  INITIAL_QUOTES=NULL, I=20L, N_POP=500,
  DIRICHLET_ALPHA=rep(1,3), RISK_PROB=c(.1,.6,.9),
  RISK_COST_PROB=c(.2,.5,.8), ALPHA=10, BETA_PRICE=1,
  V_INSURER=c(A=0,B=0,C=0), ETA_CAT=1, EPS_SHARE=.01,
  CUMULATIVE_PRESENCE=FALSE, EPS_CUMULATIVE=.25,
  LAMBDA_REALLOC=.35, DAMPING_NU=.7, DENOM_EPSILON=0,
  PRICE_LB=.01, PRICE_UB=Inf, CHARGE_FACTOR=1.02,
  T_MAX=200L, TOL_CONVERGENCE=1e-7,
  FEASIBILITY_GUARD=FALSE, GUARD_TOL=1e-8,
  GUARD_MAX_HALVINGS=26L, GUARD_MAX_ZERO_STEPS=5L,
  JOINT_REPAIR=FALSE, REPAIR_MAXEVAL=600L, REPAIR_NSTARTS=4L,
  N_STARTS=5L, quiet=TRUE
) {
  stopifnot(requireNamespace("nloptr",quietly=TRUE))
  if (!"deprecatedBehavior" %in% names(formals(nloptr::slsqp)))
    stop("Use an nloptr version supporting deprecatedBehavior=FALSE")
  scenario <- match.arg(scenario)
  names_expected <- c("A","B","C")
  if (!identical(names(ARCHETYPES),names_expected) ||
      any(!ARCHETYPES %in% c("national","communal","communal_risk")))
    stop("ARCHETYPES must be a named vector A/B/C with valid partitions")
  if (LAMBDA_REALLOC<=0 || LAMBDA_REALLOC>1 || BETA_PRICE<=0 ||
      EPS_SHARE<0 || EPS_SHARE>1 || EPS_CUMULATIVE<0 ||
      EPS_CUMULATIVE>1 || DAMPING_NU<=0 || DAMPING_NU>1)
    stop("Invalid model parameters")
  if (length(RISK_PROB)!=length(RISK_COST_PROB) ||
      length(RISK_PROB)!=length(DIRICHLET_ALPHA))
    stop("Risk-class vectors must have the same length")

  insurers <- names_expected; archetypes <- ARCHETYPES
  R <- length(RISK_PROB)
  set.seed(seed)
  N <- rep(N_POP,I)
  dirraw <- matrix(rgamma(I*R,shape=rep(DIRICHLET_ALPHA,each=I)),I,R)
  risk_props <- dirraw/rowSums(dirraw)
  u <- as.vector(risk_props %*% RISK_PROB)
  N_risk <- matrix(N,I,R) * risk_props
  risk_cost <- ETA_CAT * ALPHA * RISK_COST_PROB
  pure_ir <- outer(u,risk_cost * RISK_PROB)
  loaded_fair_ir <- CHARGE_FACTOR*pure_ir
  p0 <- CHARGE_FACTOR*sum(N_risk*pure_ir)/sum(N_risk)
  v <- as.numeric(V_INSURER[insurers]);v[!is.finite(v)]<-0
  names(v)<-insurers
  monitor <- switch(scenario,unconstrained="none",national="national",
                    commune="commune",risk_class="risk_class")

  dimension <- function(type) switch(type,national=1L,communal=I,
                                    communal_risk=I*R)
  encode <- function(mat,type) {
    mat <- matrix(as.numeric(mat),I,R)
    if (type=="national") return(mean(mat))
    if (type=="communal") return(rowMeans(mat))
    if (type=="communal_risk") return(as.vector(t(mat)))
    stop("Unknown architecture")
  }
  decode <- function(z,type) {
    if (type=="national") return(matrix(as.numeric(z)[1],I,R))
    if (type=="communal") return(matrix(rep(as.numeric(z),each=R),I,R,byrow=TRUE))
    if (type=="communal_risk") return(matrix(as.numeric(z),I,R,byrow=TRUE))
    stop("Unknown architecture")
  }
  initial <- if (is.null(INITIAL_QUOTES)) setNames(rep(p0,3),insurers)
             else INITIAL_QUOTES
  if (!identical(names(initial),insurers) || any(!is.finite(initial)) ||
      any(initial<=0)) stop("INITIAL_QUOTES must be positive named A/B/C")
  make_initial_prices <- function() setNames(lapply(insurers,function(nm)
    decode(rep(initial[[nm]],dimension(archetypes[[nm]])),archetypes[[nm]])),insurers)
  shares_target <- function(prices) {
    util <- lapply(insurers,function(nm)
      matrix(ALPHA*u,I,R)-BETA_PRICE*prices[[nm]]+v[[nm]])
    # Stable softmax; the common location utility cancels within each cell.
    umax <- Reduce(pmax,util)
    # Preserve the I x R cell structure after clipping logit utilities.
    # pmin(500, matrix) may strip dimensions because its first arg is scalar.
    numer <- lapply(util, function(z) {
      centered <- z - umax
      clipped <- pmax(pmin(centered, 500), -700)
      matrix(exp(clipped), nrow=I, ncol=R)
    })
    denom <- Reduce(`+`,numer)
    if (DENOM_EPSILON!=0) {
      # Absolute outside-option intercept in the original utility scale.
      denom <- denom + DENOM_EPSILON*exp(-umax)
    }
    result <- setNames(lapply(numer, function(z) z / denom), insurers)
    if (!all(vapply(result, function(mat)
           identical(dim(mat), c(I, R)), logical(1)))) {
      stop("shares_target(): missing I x R matrix dimensions")
    }
    result
  }
  aggregate_share <- function(mat,unit) {
    if (unit=="national") return(sum(N_risk*mat)/sum(N_risk))
    if (unit=="commune") return(rowSums(risk_props*mat))
    if (unit=="risk_class") return(colSums(N_risk*mat)/colSums(N_risk))
    stop("Invalid monitoring unit")
  }

  prices_prev <- make_initial_prices()
  shares_initial <- shares_target(prices_prev)
  shares_prev <- shares_initial
  pbar_prev <- prices_prev
  price_history <- list(prices_prev)
  share_history <- list(shares_prev)
  pbar_history <- list(pbar_prev)
  con_records <- list();opt_records<-list();guard_records<-list()
  cv_records <- list(data.frame(period=0L,max_price_change=NA_real_,
                                max_share_change=NA_real_))
  zero_steps <- 0L

  predicted_shares <- function(prices,previous_shares) {
    sig <- shares_target(prices)
    setNames(lapply(insurers,function(nm)
      (1-LAMBDA_REALLOC)*previous_shares[[nm]]+LAMBDA_REALLOC*sig[[nm]]),insurers)
  }
  share_slacks <- function(old_sh,new_sh,initial_sh,period) {
    if (monitor=="none") return(data.frame())
    rows<-list()
    for (nm in insurers) {
      old<-aggregate_share(old_sh[[nm]],monitor)
      cur<-aggregate_share(new_sh[[nm]],monitor)
      ini<-aggregate_share(initial_sh[[nm]],monitor)
      types<-if (CUMULATIVE_PRESENCE) c("period","cumulative") else "period"
      for (kind in types) {
        floor<-if(kind=="period") (1-EPS_SHARE)*old else
                (1-EPS_CUMULATIVE)*ini
        ss<-as.numeric(cur-floor)
        rows[[length(rows)+1L]]<-data.frame(period=period,insurer=nm,
          constraint_type=monitor,constraint_label=kind,unit=seq_along(ss),
          old_share=as.numeric(old),initial_share=as.numeric(ini),
          new_share=as.numeric(cur),floor=as.numeric(floor),slack=ss,
          binding=ss<1e-6)
      }
    }
    do.call(rbind,rows)
  }
  expected_profit <- function(nm,z,old_prices,old_sh,old_pbar) {
    trial<-old_prices
    trial[[nm]]<-decode(z,archetypes[[nm]])
    sig<-shares_target(trial)[[nm]]
    inherited <- (1-LAMBDA_REALLOC)*old_sh[[nm]]*(old_pbar[[nm]]-pure_ir)
    fresh <- LAMBDA_REALLOC*sig*(trial[[nm]]-pure_ir)
    sum(N_risk*(inherited+fresh))
  }
  candidate_slacks <- function(nm,z,old_prices,old_sh,period) {
    if(monitor=="none")return(0)
    trial<-old_prices;trial[[nm]]<-decode(z,archetypes[[nm]])
    newer<-predicted_shares(trial,old_sh)
    ss<-share_slacks(old_sh,newer,shares_initial,period)
    ss$slack[ss$insurer==nm]
  }
  price_box <- function(nm) {
    n<-dimension(archetypes[[nm]])
    list(lower=rep(PRICE_LB,n),upper=rep(if(is.finite(PRICE_UB))
                                        PRICE_UB else 1e6,n))
  }
  optimise_one <- function(nm,t) {
    typ <- archetypes[[nm]];bb<-price_box(nm)
    x0<-encode(prices_prev[[nm]],typ)
    starts<-list(x0,encode(loaded_fair_ir,typ))
    if (N_STARTS>2)for(k in seq_len(N_STARTS-2))
      starts[[length(starts)+1L]]<-x0*exp(rnorm(length(x0),sd=.1))
    starts<-lapply(starts,function(z)pmin(bb$upper,pmax(bb$lower,as.numeric(z))))
    objective<-function(z)-expected_profit(nm,z,prices_prev,shares_prev,pbar_prev)
    best<-NULL
    for(k in seq_along(starts)){
      arg<-list(x0=starts[[k]],fn=objective,lower=bb$lower,upper=bb$upper,
                nl.info=FALSE,deprecatedBehavior=FALSE,
                control=list(maxeval=2000,xtol_rel=1e-8))
      if(monitor!="none")
        arg$hin<-function(z) -candidate_slacks(nm,z,prices_prev,shares_prev,t)
      sol<-tryCatch(do.call(nloptr::slsqp,arg),error=function(e)NULL)
      if(is.null(sol)||is.null(sol$par)||!is.finite(sol$value)||
         is.na(sol$convergence)||sol$convergence<0)next
      if(monitor!="none"&&any(candidate_slacks(nm,sol$par,prices_prev,shares_prev,t) < -1e-6))next
      vv<- -objective(sol$par)
      if(is.null(best)||vv>best$value)
        best<-list(par=sol$par,value=vv,status=sol$convergence,start=k)
    }
    if(is.null(best)){
      warning(sprintf("SQP fallback for %s at iteration %d",nm,t))
      best<-list(par=x0,value=-objective(x0),status=NA_integer_,start=NA_integer_)
    }
    opt_records[[length(opt_records)+1L]] <<- data.frame(
      period=t,insurer=nm,n_starts=length(starts),selected_start=best$start,
      convergence=best$status,objective=best$value,
      price_min=min(best$par),price_max=max(best$par),
      hit_lower=any(abs(best$par-bb$lower)<1e-6),
      hit_upper=any(abs(best$par-bb$upper)<1e-6))
    decode(best$par,typ)
  }
  joint_vector<-function(P)unlist(lapply(insurers,function(nm)
    encode(P[[nm]],archetypes[[nm]])),use.names=FALSE)
  from_joint<-function(z){
    ans<-list();k<-0L
    for(nm in insurers){d<-dimension(archetypes[[nm]])
      ans[[nm]]<-decode(z[seq.int(k+1L,k+d)],archetypes[[nm]])
      k<-k+d
    }
    ans
  }
  evaluate_joint<-function(P,t){
    sh<-predicted_shares(P,shares_prev)
    sl<-share_slacks(shares_prev,sh,shares_initial,t)
    list(prices=P,sigma=shares_target(P),shares=sh,
         min_slack=if(nrow(sl))min(sl$slack)else Inf,
         slacks=sl)
  }
  joint_repair<-function(proposed,previous,t){
    target<-joint_vector(proposed);base<-joint_vector(previous)
    weight<-unlist(lapply(insurers,function(nm)
      switch(archetypes[[nm]],national=I*R,communal=rep(R,I),
             communal_risk=rep(1,I*R))),use.names=FALSE)
    weight<-weight/sum(weight)
    bb<-lapply(insurers,price_box)
    lb<-unlist(lapply(bb,`[[`,"lower"),use.names=FALSE)
    ub<-unlist(lapply(bb,`[[`,"upper"),use.names=FALSE)
    obj<-function(z)sum(weight*(z-target)^2)
    hfun<-function(z)-evaluate_joint(from_joint(z),t)$slacks$slack
    starts<-list(target,base,(target+base)/2,joint_vector(make_initial_prices()))
    best<-NULL
    for(j in seq_len(min(REPAIR_NSTARTS,length(starts)))){
      arg<-list(x0=pmin(ub,pmax(lb,starts[[j]])),fn=obj,hin=hfun,
                lower=lb,upper=ub,nl.info=FALSE,deprecatedBehavior=FALSE,
                control=list(maxeval=REPAIR_MAXEVAL,xtol_rel=1e-9))
      sol<-tryCatch(do.call(nloptr::slsqp,arg),error=function(e)NULL)
      if(is.null(sol)||is.null(sol$par))next
      trial<-evaluate_joint(from_joint(sol$par),t)
      if(is.finite(trial$min_slack)&&trial$min_slack>=-GUARD_TOL){
        best<-list(accepted=trial,distance=obj(sol$par),
                   adjustment=max(abs(sol$par-target)),status=sol$convergence)
        break
      }
    }
    best
  }

  for(t in seq_len(T_MAX)){
    cand_prices<-setNames(lapply(insurers,function(nm)optimise_one(nm,t)),insurers)
    candidate_for_scale<-function(scale){
      trial<-setNames(lapply(insurers,function(nm)
        prices_prev[[nm]]+DAMPING_NU*scale*(cand_prices[[nm]]-prices_prev[[nm]])),insurers)
      evaluate_joint(trial,t)
    }
    raw<-candidate_for_scale(1)
    accepted<-raw; step_scale<-1;halvings<-0L
    repaired<-FALSE;adjustment<-0;repair_distance<-NA_real_;repair_status<-NA_integer_
    if(FEASIBILITY_GUARD && raw$min_slack < -GUARD_TOL){
      accepted<-NULL
      for(k in seq_len(GUARD_MAX_HALVINGS)){
        trial<-candidate_for_scale(2^(-k))
        if(trial$min_slack >= -GUARD_TOL){
          accepted<-trial;step_scale<-2^(-k);halvings<-k;break
        }
      }
      if(is.null(accepted)&&JOINT_REPAIR){
        rescue<-joint_repair(raw$prices,prices_prev,t)
        if(!is.null(rescue)){
          accepted<-rescue$accepted;repaired<-TRUE
          step_scale<-NA_real_;adjustment<-rescue$adjustment
          repair_distance<-rescue$distance;repair_status<-rescue$status
        }
      }
      if(is.null(accepted)){
        zero<-candidate_for_scale(0)
        if(zero$min_slack >= -GUARD_TOL){
          accepted<-zero;step_scale<-0;halvings<-GUARD_MAX_HALVINGS+1L
        }else stop(sprintf("Feasible update not found at iteration %d",t))
      }
    }
    if(FEASIBILITY_GUARD && accepted$min_slack < -GUARD_TOL)
      stop("Internal feasibility guard error")
    if(is.finite(step_scale) && step_scale==0L)zero_steps<-zero_steps+1L
    else zero_steps<-0L
    if(zero_steps>=GUARD_MAX_ZERO_STEPS)stop("Safeguard stalled with repeated zero steps")
    guard_records[[length(guard_records)+1L]]<-data.frame(
      period=t,step_scale=step_scale,halvings=halvings,
      min_slack_unchecked=raw$min_slack,min_slack_accepted=accepted$min_slack,
      guarded=is.na(step_scale)||step_scale<1,
      zero_step=is.finite(step_scale)&&step_scale==0,
      joint_repair=repaired,repair_max_abs_adjustment=adjustment,
      repair_distance=repair_distance,repair_solver_status=repair_status)
    prices_new<-accepted$prices
    sigma_new<-shares_target(prices_new)
    shares_new<-setNames(lapply(insurers,function(nm)
      (1-LAMBDA_REALLOC)*shares_prev[[nm]]+LAMBDA_REALLOC*sigma_new[[nm]]),insurers)
    pbar_new<-setNames(lapply(insurers,function(nm){
      numerator<-(1-LAMBDA_REALLOC)*shares_prev[[nm]]*pbar_prev[[nm]]+
        LAMBDA_REALLOC*sigma_new[[nm]]*prices_new[[nm]]
      ifelse(shares_new[[nm]]>1e-14,numerator/shares_new[[nm]],prices_new[[nm]])
    }),insurers)
    pr_res<-max(vapply(insurers,function(nm)
      max(abs(prices_new[[nm]]-prices_prev[[nm]])),numeric(1)))
    sh_res<-max(vapply(insurers,function(nm)
      max(abs(shares_new[[nm]]-shares_prev[[nm]])),numeric(1)))
    cv_records[[length(cv_records)+1L]]<-data.frame(
      period=t,max_price_change=pr_res,max_share_change=sh_res)
    sdiag<-share_slacks(shares_prev,shares_new,shares_initial,t)
    if(nrow(sdiag))con_records[[length(con_records)+1L]]<-sdiag
    price_history[[length(price_history)+1L]]<-prices_new
    share_history[[length(share_history)+1L]]<-shares_new
    pbar_history[[length(pbar_history)+1L]]<-pbar_new
    prices_prev<-prices_new;shares_prev<-shares_new;pbar_prev<-pbar_new
    if(!quiet && (t==1||t%%10==0))message("iteration=",t,
      " quote residual=",format(pr_res,digits=4),
      " share residual=",format(sh_res,digits=4))
    if(pr_res<TOL_CONVERGENCE && sh_res<TOL_CONVERGENCE)break
  }

  # The long-form period/cell record used by the paper and the response letter.
  rows<-vector("list",length(price_history)*I*R);k<-0L
  for(tt in seq_along(price_history)){
    P<-price_history[[tt]];S<-share_history[[tt]];PB<-pbar_history[[tt]]
    for(m in seq_len(I))for(r in seq_len(R)){
      k<-k+1L
      row<-data.frame(period=tt-1L,commune_id=m-1L,risk_class=r-1L,
        commune_name=paste0("commune_",m-1L),N_commune=N[m],
        risk_proportion=risk_props[m,r],risk_prob=u[m],
        risk_cost=risk_cost[r],risk_frequency_class=RISK_PROB[r],
        segment_mass=N_risk[m,r],pure_cost=pure_ir[m,r],
        loaded_fair_premium=loaded_fair_ir[m,r],scenario_key=scenario)
      for(nm in insurers){
        exposure<-N_risk[m,r]*S[[nm]][m,r]
        row[[paste0(nm,"_prix")]]<-P[[nm]][m,r]
        row[[paste0(nm,"_pbar")]]<-PB[[nm]][m,r]
        row[[paste0(nm,"_part")]]<-S[[nm]][m,r]
        row[[paste0(nm,"_exposure")]]<-exposure
        row[[paste0(nm,"_claims")]]<-exposure*pure_ir[m,r]
        row[[paste0(nm,"_premium_revenue")]]<-exposure*P[[nm]][m,r]
        row[[paste0(nm,"_profit_contrib")]]<-exposure*(P[[nm]][m,r]-pure_ir[m,r])
        row[[paste0(nm,"_transfer_contrib")]]<-
          exposure*(P[[nm]][m,r]-loaded_fair_ir[m,r])
      }
      rows[[k]]<-row
    }
  }
  out<-do.call(rbind,rows)
  attr(out,"constraint_diagnostics")<-if(length(con_records))do.call(rbind,con_records)else data.frame()
  attr(out,"optimization_diagnostics")<-if(length(opt_records))do.call(rbind,opt_records)else data.frame()
  attr(out,"convergence_diagnostics")<-do.call(rbind,cv_records)
  attr(out,"guard_diagnostics")<-if(length(guard_records))do.call(rbind,guard_records)else data.frame()
  attr(out,"parameters")<-list(scenario=scenario,seed=seed,I=I,N_POP=N_POP,
    ARCHETYPES=ARCHETYPES,BETA_PRICE=BETA_PRICE,ETA_CAT=ETA_CAT,
    EPS_SHARE=EPS_SHARE,EPS_CUMULATIVE=EPS_CUMULATIVE,
    LAMBDA_REALLOC=LAMBDA_REALLOC,DAMPING_NU=DAMPING_NU,
    T_MAX=T_MAX,TOL_CONVERGENCE=TOL_CONVERGENCE)
  class(out)<-c("climate_granularity_revised",class(out))
  out
}

constraint_diagnostics<-function(x)attr(x,"constraint_diagnostics")
optimization_diagnostics<-function(x)attr(x,"optimization_diagnostics")
guard_diagnostics<-function(x)attr(x,"guard_diagnostics")
convergence_diagnostics<-function(x)attr(x,"convergence_diagnostics")
parameter_table<-function(x)attr(x,"parameters")

# Check invariants for every accepted trajectory, not just selected endpoints.
check_market_invariants<-function(x,tol=1e-8){
  stopifnot(inherits(x,"climate_granularity_revised"))
  nm<-c("A","B","C")
  for(i in nm)stopifnot(all(is.finite(x[[paste0(i,"_part")]])))
  s<-Reduce(`+`,lapply(nm,function(i)x[[paste0(i,"_part")]]))
  if(max(abs(s-1))>tol)stop("Shares do not sum to one")
  sl<-constraint_diagnostics(x)
  if(nrow(sl) && min(sl$slack) < -tol)stop("Reported trajectory violates monitoring constraints")
  invisible(TRUE)
}

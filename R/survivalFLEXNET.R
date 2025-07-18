survivalFLEXNET <- function(formula, data, ratetable, m=3, mpos = NULL, mquant = NULL, init = NULL, 
                             delta_th = 0, weights=NULL)
{
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a ratetable argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  if(!is.null(init)){if(!is.numeric(init))stop("Argument 'init' must be a vector of numeric values.") }
  if(!is.numeric(delta_th))stop("'delta_th' must be numeric.")
  if(length(delta_th) != 1) stop("'delta_th' must be a single value.") 
  if(!is.null(mpos) & !is.null(mquant))warning("'mpos' and 'mquant' have both been specified. 'mpos' values have been chosen over 'mquant' quantiles.") #ligne : mquant = NULL
  
  ####### data management
  
  time <- data[,as.character(formula[[2]][2])] 
  event <- data[,as.character(formula[[2]][3])]
  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms)>1) stop("More than one 'strata' term found in  the formula. Only one variable at a time can be stratified")
  ratetable_terms <- grep("^ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  if(length(ratetable_terms)>1) stop("More than one 'ratetable' term found in  the formula.")
  CV <- setdiff(all_terms, c(strata_terms, ratetable_terms))
  if(length(CV) == 0){covnames = "1"} else{covnames <- CV}  
  label<- NULL
  covs <- as.formula(paste("~", paste(covnames, collapse = " + ")))
  cova <- model.matrix(covs, data)[, -1, drop = FALSE]
  covnames <- colnames(cova)
  extract_vars <- function(term) {
    var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
    vars <- trimws(unlist(strsplit(var_string, ",")))
    return(vars)
  }
  assign_ratetable_vars <- function(vars) {
    age <- year <- sex <- NULL
    for (var in vars) {
      if (grepl("age = ", var)) {
        age <- sub("age = ", "", var)
      } else if (grepl("year = ", var)) {
        year <- sub("year = ", "", var)
      } else if (grepl("sex = ", var)) {
        sex <- sub("sex = ", "", var)
      }
    }
    unnamed_vars <- setdiff(vars, c(age, sex, year))
    if (length(unnamed_vars) > 0) {
      if (is.null(age) && length(unnamed_vars) >= 1) age <- unnamed_vars[1]
      if (is.null(year) && length(unnamed_vars) >= 2) year <- unnamed_vars[2]
      if (is.null(sex) && length(unnamed_vars) >= 3) sex <- unnamed_vars[3]
    }
    return(list(age = age, year = year, sex = sex))
  }
  strata_var = unlist(lapply(strata_terms, extract_vars))
  if(!is.null(strata_var) && strata_var %in% covnames) stop("The stratified covariate also appears as a covariate in the formula.")
  if(is.null(strata_var)){
    timevar = strata_var
    xlevels = NULL
  }
  if(!is.null(strata_var)){
    timevar <- data[,strata_var]
    xlevels <- list(levels(as.factor(timevar)))
    names(xlevels) <- c(strata_var)
  }
  if(!is.null(timevar)){
    if(length(unique(timevar))>15) stop("The variable with a time-dependant effect has too many categories (>15)")
  }
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
  age <- data[,ratetable_vars$age]  
  year <- data[,ratetable_vars$year]
  sex <- data[,ratetable_vars$sex] 
  if(is.null(weights)){weights = rep(1,dim(data)[1])}
  if(max(time)<30) stop("The event times must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.numeric(age)) stop("'age' must be numeric")
  if(max(age)<30) stop("The ages must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.character(sex)) stop("'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("'sex' must be 'male' or 'female'")
  # if(!is.date(year)) stop("The values for 'year' must be of the 'date' class")
  if(!(is.numeric(year) || inherits(year, "Date") || inherits(year, "POSIXct") || inherits(year, "POSIXlt")))stop("'year' must be of class numeric, Date, POSIXct, or POSIXlt.")
  if(!is.null(weights)){
    if(!is.numeric(weights)) stop("Argument 'weights' must be a numeric vector")
    if(length(weights)!=dim(data)[1]) stop("Argument 'weights' must have the same length as the number of rows of the 'data' argument. (", dim(data)[1],")")}
  if(!is.null(mpos) & !is.null(mquant)){mquant = NULL}
  d <- cbind(time, event, cova, age, 1*(sex=="male"), year)
  na <- !is.na(apply(d, MARGIN=1, FUN = "sum"))
  
  time <- time[na]
  event <- event[na]
  cova <- as.matrix(cova[na,])
  age <- age[na]
  sex <- sex[na]
  year <- year[na]
  
  ###### Compute the expected mortality of individuals
  
  
  hP <- sapply(seq_along(time), function(i) {
    expectedhaz(ratetable = ratetable, age = age[i], sex = sex[i], year = year[i], time = time[i])
  })
  
  ###### log likelihood functions
  
  if (!is.null(covnames) & is.null(timevar)){
    
    logll1 <- function(beta, gamma, time, event, cova, hP, w, m, mpos, mquant){
      return(-1*sum( w*(
        event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos, mquant)$spln *
                      exp(splinecube(time, gamma, m, mpos, mquant)$spln + cova %*% beta) ) -
          exp(splinecube(time, gamma, m, mpos, mquant)$spln + cova %*% beta)
      ) ) ) } ###p5489 StatinMed P.Nelson PC. Lambert Flexible Parametric models for relative survival
    
    
    if(m == 0){gamma_names = NULL
    }else{gamma_names <- paste0("gamma", 2:(m+1))}
    
    label <- c(covnames, "gamma0",
               "gamma1", gamma_names)
    
    if(!is.null(init)){
      if(length(init) != dim(cova)[2]+m+2) stop("'init' length must be ", dim(cova)[2]+m+2,
                                                " ( ", dim(cova)[2]," covariate(s) and ",m+2," parameters for the Restricted Cubic Spline).")
      init1 <- init
    }else{init1 <- c(rep(0,dim(cova)[2]+m+2))}
    
    loglik1 <- function(par, time, event, cova, hP, w, m, mpos, mquant){
      beta <- par[1:dim(cova)[2]]
      gamma <- par[(dim(cova)[2]+1):length(par)]
      return(logll1(beta, gamma, time, event, cova, hP, w, m, mpos, mquant)) }
    
    suppressWarnings({
      logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                         cova = cova, hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
      
      indic <- 0
      while(indic <= 5){
        ll_val <- logllmax1$value
        logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                           event = event, cova = cova, hP = hP, w = weights,
                           m = m, mpos = mpos, mquant = mquant)
        delta <- ll_val - logllmax1$value
        if(delta_th == 0){
          if(delta == delta_th) {indic = indic + 1}
        }else{ 
          if(0 < delta & delta <= delta_th) {indic = indic + 1}
        }      }
      
      logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                         event = event, cova = cova, hP = hP, w = weights,
                         m = m, mpos = mpos, mquant = mquant,
                         hessian = TRUE)
    })
    
  }
  
  if(!is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    K <- sort(unique(timevarnum))
    
    logll2 <- function(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, mquant, k) {
      
      betak <- beta
      gammak <- gamma
      idx <- covatime == k
      
      timek <- time[idx]
      eventk <- event[idx]
      hPk <- hP[idx]
      covak <- cova[idx, , drop = FALSE]
      wk <- w[idx]
      
      splk <- splinecube(timek, gammak, m, mpos, mquant)$spln
      splkP <- splinecubeP(timek, gammak, m, mpos, mquant)$spln
      linpred <- splk + as.matrix(covak) %*% as.matrix(betak)
      
      value_strate <- -1 * sum(wk * (eventk * log(hPk + (1 / timek) * splkP * exp(linpred)) - exp(linpred)))
      
      
      return(value_strate)
    }
    
    label <- c()
    for (k in seq_along(K)) {
      label <- c(label,
                 paste0(covnames, "_", k),
                 paste0("gamma", K[k], "_", 0:(m + 1)))
    }
    
    # Parameter initialization
    n_beta <- length(K) *dim(cova)[2]
    n_gammak <- length(K) * (m + 2)
    
    if (!is.null(init)) {
      if (length(init) != n_beta + n_gammak)
        stop("'init' length must be ", n_beta + n_gammak, " (",
             n_beta, " covariate(s), ", dim(cova)[2]," per level of the time-dependant covariate and ", n_gammak, " splines parameters, ", m+2," per level of the time-dependant covariate).")
      init1 <- init
    } else {
      init1 <- rep(0, n_beta +  n_gammak)
    }
    
    loglik2 <- function(par, time, event, cova, covatime, hP, w, m, mpos, mquant, kn) {
      beta <- par[1:dim(cova)[2]]
      gamma <- par[(dim(cova)[2] + 1):length(par)]
      
      logll2(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, mquant, kn)
    }
    
    
    init1 <- matrix(init1, nrow = dim(cova)[2]+m+2)
    results_ll <-list()
    for(k in K){
      
      suppressWarnings({
        logllmax1 <- optim(par = init1[,k], fn = loglik2, time = time, event = event,
                           cova = cova, covatime = timevarnum, hP = hP, w = weights,
                           m = m, mpos = mpos, mquant = mquant, kn = K[k])
        
        indic <- 0
        while (indic <= 5) {
          last_logllmax1 <- logllmax1
          ll_val <- logllmax1$value
          
          logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, event = event,
                             cova = cova, covatime = timevarnum, hP = hP, w = weights,
                             m = m, mpos = mpos, mquant = mquant, kn = K[k])
          
          delta <- ll_val - logllmax1$value
          if (delta_th == 0) {
            if (delta == delta_th) indic <- indic + 1
          } else {
            if (0 < delta & delta <= delta_th) indic <- indic + 1
          }
        }
        
        tryCatch({
          logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, event = event,
                             cova = cova, covatime = timevarnum, hP = hP, w = weights,
                             m = m, mpos = mpos, mquant = mquant, kn = K[k],
                             hessian = TRUE)
        }, error = function(e) {
          logllmax1 <- last_logllmax1
        })
      })
      results_ll[[k]] <- logllmax1
    }
  }
  
  if (is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    K = sort(unique(timevarnum))
    
    logll3 <- function(gamma, time, event, covatime, hP, w, m, mpos, mquant, k){
      
      gammak <- gamma
      idx <- covatime == k
      
      timek <- time[idx]
      eventk <- event[idx]
      hPk <- hP[idx]
      covak <- cova[idx, , drop = FALSE]
      wk <- w[idx]
      
      splk <- splinecube(timek, gammak, m, mpos, mquant)$spln
      splkP <- splinecubeP(timek, gammak, m, mpos, mquant)$spln
      linpred <- splk 
      
      value_strate <- -1*sum(wk * (eventk * log(hPk + (1/timek)* splkP *
                                                  exp(linpred ) ) - exp(linpred)
      )
      )
      return(value_strate)
      
    }
    
    
    label <- c()
    for (k in seq_along(K)) {
      label <- c(label, paste0("gamma", K[k], "_", 0:(m + 1)))
    }
    
    
    # Parameter initialization
    n_gammak <- length(K) * (m + 2)
    
    if (!is.null(init)) {
      if (length(init) != n_gammak)
        stop("'init' length must be ",  n_gammak, " (",  n_gammak, " spline parameters per level of the time-dependent variable).")
      init1 <- init
    } else {
      init1 <- rep(0, n_gammak)
    }
    
    loglik3 <- function(par, time, event, covatime, hP, w, m, mpos, mquant, kn){
      
      gamma <- par
      
      return(logll3(gamma, time, event, covatime, hP, w, m, mpos, mquant, kn)) }
    
    init1 <- matrix(init1, nrow = m+2)
    results_ll <-list()
    
    for(k in K){
      
      suppressWarnings({
        logllmax1 <- optim(par = init1[,k], fn = loglik3, time = time, event = event,
                           covatime = timevarnum, hP = hP, w = weights
                           , m = m, mpos = mpos, mquant = mquant, kn = K[k])
        
        indic <- 0
        while(indic <= 5){
          last_logllmax1 <- logllmax1
          ll_val <- logllmax1$value
          
          logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                             event = event, covatime = timevarnum,
                             hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, kn = K[k])
          delta <- ll_val - logllmax1$value
          if(delta_th == 0){
            if(delta == delta_th) {indic = indic + 1}
          }else{ 
            if(0 < delta & delta <= delta_th) {indic = indic + 1}
          }        }
        
        tryCatch({logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                                     event = event, covatime = timevarnum,
                                     hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, kn = K[k],
                                     hessian = TRUE)}
                 , error = function(e){logllmax1 <- last_logllmax1})
      })
      
      results_ll[[k]] <- logllmax1
      
    }
    ## end K loop
  }
  
  #NULL model
  logll0 <- function(gamma, time, event, hP, w, m, mpos, mquant){
    return(-1*sum(w*(
      event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos, mquant)$spln *
                    exp(splinecube(time, gamma, m, mpos, mquant)$spln) ) -
        exp(splinecube(time, gamma, m, mpos, mquant)$spln)
    ) )
    ) 
  }
  
  if(m == 0){gamma_names = NULL
  }else{gamma_names <- paste0("gamma", 2:(m+1))}
  
  labelNULL <- c(covnames, "gamma0",
                 "gamma1", gamma_names)
  
  if(is.null(covnames) & is.null(timevar)){
    if(!is.null(init)){
      if(length(init) != (m+2)) stop("'init' length must be ",(m+2)," (",(m+2)," parameters for the Restricted Cubic Spline).")
      init0 <- init}else{init0 <- rep(0,m+2)
      }}else{
        init0 <- rep(0,m+2)
      }
  
  
  loglik0 <- function(par, time, event, hP, w, m, mpos, mquant){
    gamma <- par
    return(logll0(gamma, time, event, hP, w, m, mpos, mquant)) }
  
  suppressWarnings({
    
    logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                       hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
    
    indic <- 0
    while(indic <= 5){
      ll_val <- logllmax0$value
      logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                         event = event, hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
      delta <- ll_val - logllmax0$value
      if(delta_th == 0){
        if(delta == delta_th) {indic = indic + 1}
      }else{ 
        if(0 < delta & delta <= delta_th) {indic = indic + 1}
      }  
    }
    
    logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                       event = event, hP = hP, w = weights, 
                       hessian = TRUE, m = m, mpos = mpos, mquant = mquant)
  })
  
  ##récupération de mpos et mquant si ils n'ont pas été spécifiés dans la formule.
  
  if(is.null(timevar)){
    if(is.null(mpos)){
      if(is.null(mquant)){
        a <- c()
        for(i in (0:(m+1))){
          a <- c(a,i/(m+1))}
        mpos <- quantile(log(time), probs = a)
        mpos <- as.numeric(mpos)
        mquant <- a 
      }else{
        a <- c(mquant)
        mpos <- quantile(log(time), probs = a)
      }
    }
  }
  if(!(is.null(timevar))){
    if(is.null(mpos)){
      if(is.null(mquant)){
        a <- c()
        for(i in (0:(m+1))){
          a <- c(a,i/(m+1))}
        mpos_strates <- c()
        for(k in seq_along(K)){
          idx <- timevarnum == k
          timek <- time[idx]
          
          mpos_strates <- c(mpos_strates, quantile(log(timek), probs = a) )
          mquant <- a 
        }
        mpos_strates <- as.numeric(mpos_strates)
        mpos <- matrix(mpos_strates, nrow = length(K))
        
      }else{
        a <- c(mquant)
        mpos_strates <- c()
        for(k in seq_along(K)){
          idx <- timevarnum == k
          timek <- time[idx]
          
          mpos_strates <- c(mpos_strates, quantile(log(timek), probs = a) )
          mquant <- a 
        }
        mpos_strates <- as.numeric(mpos_strates)
        mpos<- matrix(mpos_strates, nrow = length(K))
      }
    }
  }
  if (!is.null(covnames) & is.null(timevar)){
    
    t.table <- data.frame(coef = logllmax1$par,
                          ecoef = exp(logllmax1$par),
                          se = sqrt(diag(solve(logllmax1$hessian))),
                          z = logllmax1$par/sqrt(diag(solve(logllmax1$hessian))),
                          p = 2*(1-pnorm(abs(logllmax1$par/sqrt(diag(solve(logllmax1$hessian)))), 0, 1)),
                          row.names = label)
  }
  
  if (!is.null(covnames) & !is.null(timevar)){
    
    ret_par <- c()
    sq_solve <- c()
    for(i in 1:length(K)){
      ret_par <- c(ret_par, results_ll[[i]]$par)
      sq_solve <- c(sq_solve, sqrt(diag(solve(results_ll[[i]]$hessian) )))
    }
    
    t.table <- data.frame(coef = ret_par,
                          ecoef = exp(ret_par),
                          se = sq_solve,
                          z = ret_par/sq_solve,
                          p = 2*(1-pnorm(abs(ret_par/sq_solve), 0, 1)),
                          row.names = label)
  }
  
  if(is.null(covnames) & is.null(timevar)){
    t.table <- data.frame(coef = logllmax0$par,
                          ecoef = exp(logllmax0$par),
                          se = sqrt(diag(solve(logllmax0$hessian))),
                          z = logllmax0$par/sqrt(diag(solve(logllmax0$hessian))),
                          p = 2*(1-pnorm(abs(logllmax0$par/sqrt(diag(solve(logllmax0$hessian)))), 0, 1)),
                          row.names = labelNULL)
  }
  
  if(is.null(covnames) & !is.null(timevar)){
    
    ret_par <- c()
    sq_solve <- c()
    for(i in 1:length(K)){
      ret_par <- c(ret_par, results_ll[[i]]$par)
      sq_solve <- c(sq_solve, sqrt(diag(solve(results_ll[[i]]$hessian) )))
    }
    
    t.table <- data.frame(coef = ret_par,
                          ecoef = exp(ret_par),
                          se = sq_solve,
                          z = ret_par/sq_solve,
                          p = 2*(1-pnorm(abs(ret_par/sq_solve), 0, 1)),
                          row.names = label)
  }
  
  names(t.table) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  
  coefficients <- t.table$coef
  names(coefficients) <- label
  
  if(is.null(timevar)){
    betaestim <- coefficients[covnames]
    lp <- cova %*% betaestim
  }
  
  if(!is.null(timevar)){
    if(!is.null(covnames)){
      pattern <- paste(covnames, collapse = "|")  
      betaestim <- coefficients[grep(pattern, names(coefficients))]
      
      lp <- numeric(nrow(cova))
      
      for (k in K) {
        idx <- which(timevarnum == k)
        
        coefnames_k <- paste0(covnames, "_", k)
        
        betak <- betaestim[coefnames_k]
        
        covak <- as.matrix(cova[idx, covnames, drop = FALSE])
        
        lp[idx] <- covak %*% betak
      }
      
    }
    if(is.null(covnames)){
      lp <- NA
    }
  }
  
  dimnames(cova)[[2]] <- covnames
  
  var <- if(!is.null(covnames) & is.null(timevar)){solve(logllmax1$hessian)
  } else if(!is.null(covnames) & !is.null(timevar)){
    res_solve <- list()
    
    for(i in 1:length(K)){
      res_solve[[i]] <- solve(results_ll[[i]]$hessian) 
    }
    res_solve
  }else if(is.null(covnames) & !is.null(timevar)) { 
    
    res_solve <- list()
    
    for(i in 1:length(K)){
      res_solve[[i]] <- solve(results_ll[[i]]$hessian) 
    }
    res_solve
    
  }else{solve(logllmax0$hessian)}
  
  
  loglik <- if(!is.null(covnames) & is.null(timevar) ){c(-1*logllmax1$value, -1*logllmax0$value)
  }else if(!is.null(covnames) & !is.null(timevar)){ 
    res_log <- c()
    
    for(i in 1:length(K)){
      res_log <- c(res_log, results_ll[[i]]$value) 
    }
    res_log <- c(-1*sum(res_log), -1*logllmax0$value, -1*res_log)
    
  }else if(is.null(covnames) & !is.null(timevar)){  
    res_log <- c()
    
    for(i in 1:length(K)){
      res_log <- c(res_log, results_ll[[i]]$value) 
    }
    res_log <- c(-1*sum(res_log), -1*logllmax0$value, -1*res_log)
  }else{-1*logllmax0$value}
  
  if(!is.null(covnames) & is.null(timevar)){names(loglik) <- c("Model", "Null model")
  }else if(!is.null(covnames) & !is.null(timevar)){ names(loglik) <- c("Model", "Null Model", paste0("Model_", timevarnames))
  
  }else if(is.null(covnames) & !is.null(timevar)){ names(loglik) <- c("Model", "Null Model", paste0("Model_", timevarnames))
  
  }else{
    names(loglik) <- c("Null Model")
  }
  
  res <- list(
    formula = formula,
    coefficients =  coefficients,
    var = var,
    t.table = t.table,
    loglik = loglik,
    linear.predictors = as.vector(lp),
    missing = !na,
    n = length(time),
    nevent = sum(event),
    y = cbind(time = time, status = event),
    x = cova,
    ays = data.frame(age = age, year = year, sex = sex),
    m = m,
    mpos = mpos,
    mquant = mquant
  )
  if (!is.null(xlevels)) {
    res$correstab <- correstab
    res$xlevels <- xlevels
    res$levelsval <- data[,strata_var]
  }
  class(res) <- "survivalNET"
  return(res)
}



dmn.em <- function(kmeans.res, Wx, bin.width, alpha, M, K, Lx,  N, verbose, 
                   maxNumOptIter, binned.data, eta, nu, etah, nuh, numOptRelTol, 
                   EM.maxit, EM.threshold, optim.options=NULL, hessian=FALSE) {
  
  
  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
  
  #lambda is log(alpha)
  alpha <- lapply(alpha, function(a){a[a <= 0] <- 1e-6;a})
  lambda <- lapply(alpha, log) #natural logarithm by default
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  
  weights <- rowSums(Ez) #initial values for \bm{\pi}
  EM.diagnostics <- data.frame()
  if (verbose) {
    cat('Expectation Maximization setup\n')
    
    pb <- txtProgressBar(min = 1,
                         max = ifelse(M*K > 1, M*K, 2),
                         initial = 1,
                         title='Numerical optimization',
                         width=40,
                         style=3)
  }
    
  #Initialize the lambda values by optimizing Q wrt the lambda values
  for (m in 1:M) { #over data types
  
    for (k in 1:K) { #over clusters
      if (verbose)
        setTxtProgressBar(pb, (m-1)*K+k)
      hkm <- vector(mode = 'list', maxNumOptIter+1) #comes from gradient computation
      gradient<- vector(mode = 'list', maxNumOptIter+1)

      #the lambda_{kj}^{(m)} can be optimized for each k and m separately
      
     optim.result <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,],
                                          hkm=hkm,
                                          gradient=gradient,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol, 
                                          hessian=hessian, 
                                          optim.options=optim.options)
      
     
     hkm <- unlist(hkm)
     gradient <- unlist(gradient)
     
    hkm <- data.frame(Datatype=names(binned.data)[m],
                         Component=k,
                         EM.iter=0,
                         hkm=hkm,
                         nll=optim.result$value,
                         gradient=gradient)
    
     if(hessian==TRUE){
       hkm$detH=det(optim.result$hessian)
     }
     
     EM.diagnostics <- rbind(EM.diagnostics, hkm)
   
    lambda[[m]][k,] <-optim.result$par
          
    }
  } #initialization ends
  
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  if (verbose)
    cat('\nExpectation Maximization\n')
  
  #EM loop
  iter <- 0
  last.nll <- .Machine$double.xmax
 
  
  real.change<- .Machine$double.xmax #1.797693e+308
    
  #use negative LL/unnormalized posterior to check the convergence
  nll_list<- vector(mode = 'list', EM.maxit+1)
  #function neg_log_likelihood is used to check the convergence
  nll_iter=1
  nll_list[[nll_iter]] <- neg_log_likelihood(weights, lambda, binned.data, eta, nu, etah, nuh)
  last.nll <- nll_list[[nll_iter]]
  nll_iter=nll_iter+1
  
  
    
  #E-step
  #M-step
  
  #while (iter < EM.maxit ) {
  while ((iter < EM.maxit) && (real.change > EM.threshold)) {
    if (verbose)
      cat('Calculating Ez values...\n')
    
    #Computes the E-step, the posterior probabilities of the cluster labels
    #The old Ez values not really used
    # returns K x N matrix Z
    Ez <- calc_z(Ez, binned.data, weights, lambda)
    stopifnot(!is.na(Ez), !is.infinite(Ez))
    
    #Computes M-step, same for lambda as the EM initialization
    if (verbose) {
      cat('Optimizing lambda...\n')
      
      pb <- txtProgressBar(min = 1,
                           max = ifelse(M*K > 1, M*K, 2),
                           initial = 1,
                           title='Numerical optimization',
                           width=40,
                           style=3)
    }
    
    for (m in seq_len(M)) {
       for (k in seq_len(K)) {
         if (verbose)
          setTxtProgressBar(pb, (m-1)*K+k)
        
        hkm <- vector(mode = 'list', maxNumOptIter+1) #comes from gradient computation
        gradient <- vector(mode = 'list', maxNumOptIter+1)
        optim.result <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,],
                                          hkm=hkm,
                                          gradient=gradient,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol, 
                                          hessian=hessian, 
                                          optim.options=optim.options)
        
        
        lambda[[m]][k,] <-optim.result$par
        hkm <- unlist(hkm)
        gradient <- unlist(gradient)
       
    
        hkm <- data.frame(Datatype=names(binned.data)[m],
                            Component=k,
                            EM.iter=iter+1,
                            hkm=hkm,
                            nll=optim.result$value,
                            gradient=gradient)
        
        
        
        if(hessian==TRUE){
          hkm$detH=det(optim.result$hessian)
        }
        
        EM.diagnostics <- rbind(EM.diagnostics, hkm)
      }
    }
    
    stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
    
        
    weights <- rowSums(Ez) # \pi_k values Ez is K x N matrix, are these normalized? Should this be rowMeans?
    
    if (verbose)
      cat('\nCalculating negative log likelihood...\n')
    
    #calculates neg. unnormalized posterior for convergence
    # weights 1xK (unnormalized, are normalized by this function/div by N)
    #lambda list of K times S matrices
    #binned.data N times S
    nll <- neg_log_likelihood(weights, lambda, binned.data, eta, nu, etah, nuh)
    stopifnot(!is.na(nll), !is.infinite(nll))
    nll_list[[nll_iter]]=nll
    nll_iter=nll_iter+1
  
    real.change=last.nll-nll
    
    if(real.change< -numOptRelTol){
      print("Warning: Neg.LL does not decrease!!!")
    }
    
    last.nll <- nll
    iter <- iter+1
    
    if (verbose){
      print(paste('--> EM Iteration:', iter, 'Neg.LL change (real):', round(real.change, -log10(numOptRelTol)) ))
      print(paste0("Neg.LL ",round(nll,-log10(numOptRelTol)) ))
    }
    
  } #EM loop ends
  
  
  # Model selection
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  
  P <- K*sum(Lx)+K-1 #k=S x K x M + (K-1) This should change for different every M_k?
  gof.BIC <- last.nll + 0.5 * log(N) * P #this is -BIC
  gof.AIC <- last.nll + P #this is 0.5*AIC
  gof <- c(NLE=last.nll, BIC=gof.BIC, AIC=gof.AIC)  #goodness of fit
  
  result <- list()
  
  result$GoodnessOfFit <- gof
  result$Ez <- Ez
  result$Group <- t(Ez)
  result$nll.data=nll.data
  result$Mixture <- list(Weight=weights/N)
  

  EM.diagnostics <- plyr::ddply(EM.diagnostics, c('Datatype', 'Component', 'EM.iter'),
                                transform, NO.iter.count=length(hkm), NO.iter=seq_along(hkm))
  result$EM.diagnostics <- EM.diagnostics
  
  if(verbose) {
    print('Mixture weights: ')
    print(weights)
    print('Hard labels:')
    print(table(apply(result$Group, 1, which.max)))
  }
  
    
  result$Fit <- list(Estimate=lapply(lambda, function(x)t(exp(x)))) #alpha parameters
  result$Data <- binned.data 
  
  result
}
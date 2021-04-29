dmn.em.flip <- function(kmeans.res, Wx, bin.width, S, alpha, M, K, Lx,  N, verbose, 
                   maxNumOptIter, binned.data, eta, nu, etah, nuh, numOptRelTol,  
                   EM.maxit, EM.threshold) {
  
  
  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
  if(shift==TRUE && flip==TRUE){
    
    La=Wx/bin.width+S-1
    #extend the length of alpha to La
    extend=floor(S/2)
    alpha_extend<-alpha
    for(m in 1:M){
      tmp=matrix(0, nrow=K, ncol=extend*2+Lx)
      for(k in 1:K){
        tmp[k,]=c(rep(0, extend), alpha[[m]][k,], rep(0, extend))
        tmp[k,1:extend]=alpha[[m]][k,1]
        tmp[k,ncol(tmp)-((extend-1):0)]=alpha[[m]][k,ncol(alpha[[m]])]
      }
      alpha_extend[[m]]=tmp/rowSums(tmp) #normalize
      
    }
    
    alpha=alpha_extend
    
    Ez2 <- array(0,dim=c(K,S,2,N))#This is 4D for each i
    for(i in 1:N){
      for(k in 1:K){
        Ez2[k,,,i]=Ez[k,i]
      }
    }
    Ez=Ez2
  }else if(shift==TRUE && flip==FALSE){ #only shift
    
  }else if(shift==FALSE && flip==TRUE){ #only flip
    
  }else{ # shift==FALSE && flip==FALSE
    
  }
  
  #lambda is log(alpha)
  alpha <- lapply(alpha, function(a){a[a <= 0] <- 1e-6;a})
  lambda <- lapply(alpha, log) #natural logarithm by default
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  
  weights <- rowSums(Ez) #initial values for \bm{\pi}
  
  if (verbose) {
    cat('Expectation Maximization setup\n')
    
    pb <- txtProgressBar(min = 1,
                         max = ifelse(M*K > 1, M*K, 2),
                         initial = 1,
                         title='Numerical optimization',
                         width=40,
                         style=3)
  }
  lambda_optim_message<-vector(mode = 'list', M)
  hkm_list<-vector(mode = 'list', M)
  
  #Initialize the lambda values by optimizing Q wrt the lambda values
  for (m in 1:M) { #over data types
    lambda_optim_message[[m]]<-list()
    hkm_list[[m]]<-list()
    for (k in 1:K) { #over clusters
      if (verbose)
        setTxtProgressBar(pb, (m-1)*K+k)
      hkm <- vector(mode = 'list', maxNumOptIter+1)
      
      #the lambda_{kj}^{(m)} can be optimized for each k and m separately
      if(shift==TRUE && flip==TRUE){
        optimise_lambda_k_shift_flip(LambdaK=lambda[[m]][k,],
                                     data=binned.data[[m]],
                                     Z=Ez[k,,,],
                                     hkm=hkm,
                                     eta=eta,
                                     nu=nu,
                                     etah=etah,
                                     nuh=nuh,
                                     verbose=verbose,
                                     MAX_GRAD_ITER=maxNumOptIter,
                                     reltol=numOptRelTol)
        
      }else if(shift==TRUE && flip==FALSE){
        
      }else if(FALSE==TRUE && flip==TRUE){
        
      }
      else{#shift==FALSE && flip==FALSE
        optim.result <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,],
                                          hkm=hkm,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol)
      }
      
      
      hkm_list[[m]][[k]]=hkm
      lambda[[m]][k,] <-optim.result$par
      lambda_optim_message[[m]][[k]]<-optim.result
      
    }
  } #initialization ends
  
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  if (verbose)
    cat('\nExpectation Maximization\n')
  
  #EM loop
  iter <- 0
  last.nll <- 0
  nll.change <- .Machine$double.xmax #1.797693e+308
  EM.diagnostics <- data.frame()
  
  EM_lambda_optim_message <- vector(mode = 'list', EM.maxit+1)
  EM_lambda_optim_message[[1]]<- lambda_optim_message
  EM_hkm_list<- vector(mode = 'list', EM.maxit+1)
  EM_hkm_list[[1]]=hkm_list
  
  nll_list<- vector(mode = 'list', EM.maxit+1)
  #function neg_log_likelihood is used to check the convergence
  nll_list[[1]] <- neg_log_likelihood(weights, lambda, binned.data, eta, nu, etah, nuh)
  nll.change_list<- vector(mode = 'list', EM.maxit+1)
  
  
  #E-step
  #M-step
  #shifting and flipping, one needs to write separate EM-loops for (shift=false, flip=false), 
  #(shift=true, flip=false), (shift=false, flip=true), (shift=true, flip=true)
  #4 different em-functions
  while ((iter < EM.maxit) && (nll.change > EM.threshold)) {
    
    if (verbose)
      cat('Calculating Ez values...\n')
    
    #Computes the E-step, the posterior probabilities of the cluster labels
    #The old Ez values not really used
    # returns K x N matrix Z
    #why offset subtracted?
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
    lambda_optim_message<-vector(mode = 'list', M)
    hkm_list<-vector(mode = 'list', M)
    for (m in seq_len(M)) {
      lambda_optim_message[[m]]<-list()
      hkm_list[[m]]<-list()
      for (k in seq_len(K)) {
        if (verbose)
          setTxtProgressBar(pb, (m-1)*K+k)
        
        hkm <- vector(mode = 'list', maxNumOptIter+1)
        optim.result <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,],
                                          hkm=hkm,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol)
        
        
        lambda[[m]][k,] <-optim.result$par
        lambda_optim_message[[m]][[k]]<-optim.result
        
        
        hkm_list[[m]][[k]] <- hkm
        
        
        hkm <- unlist(hkm)
        hkm <- data.frame(Datatype=names(binned.data)[m],
                          Component=k,
                          EM.iter=iter,
                          hkm=hkm, nll=optim.result$value)
        EM.diagnostics <- rbind(EM.diagnostics, hkm)
      }
    }
    #save.image("/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/DMM-private-master-devel-works/shif.flip.debugging.RData")
    EM_lambda_optim_message[[iter+2]]<- lambda_optim_message
    
    EM_hkm_list[[iter+2]]=hkm_list
    stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
    
    # if ((shift.ratio > 0 &&shift) || flip){
    #   if (verbose)
    #     cat('\nOptimizing shift/flip parameters...\n')
    #   #this is in util.R
    #   #optimization_func is in dmn.cpp, what does it try to optimize?
    #   #returns $shift (length N vector)
    #   #returns $flip (length N locigal)
    #   best.shift.and.flip.params <- brute.force.integer.minimizer(fn=optimization_func,
    #                                                  limits=attr(binned.data, 'shift.limits'),
    #                                                  step=bin.width,
    #                                                  parallel.shift=parallel.shift,
    #                                                  verbose=verbose,
    #                                                 shift=shift.reads,
    #                                                  flip=flip,
    #                                                  data=binned.data,
    #                                                  alpha=lapply(lambda, exp),
    #                                                  Z=Ez)
    #   if (verbose) {
    #     cat('\nRead shifting frequencies: \n')
    #     print(table(best.shift.and.flip.params$shift))
    # 
    #     cat('\nRead flip parameters: \n')
    #     print(table(best.shift.and.flip.params$flip))
    # 
    #   }
    #   binned.data <- shift.and.flip.signal(binned.data, seq_len(N),
    #                                        best.shift.and.flip.params$shift,
    #                                        best.shift.and.flip.params$flip)
    # } #end of shifting and flipping
    
    weights <- rowSums(Ez) # \pi_k values Ez is K x N matrix, are these normalized? Should this be rowMeans?
    
    if (verbose)
      cat('\nCalculating negative log likelihood...\n')
    
    #calculates neg. unnormalized posterior for convergence
    # weights 1xK (unnormalized, are normalized by this function/div by N)
    #lambda list of K times S matrices
    #binned.data N times S
    nll <- neg_log_likelihood(weights, lambda, binned.data, eta, nu, etah, nuh)
    stopifnot(!is.na(nll), !is.infinite(nll))
    nll_list[[iter+2]]<-nll
    
    nll.change <- abs(last.nll - nll)
    nll.change_list[[iter+1]]<-nll.change
    last.nll <- nll
    
    iter <- iter+1
    
    if (verbose)
      print(paste('--> EM Iteration:', iter, 'Neg.LL change:', round(nll.change, 6)))
  } #EM loop ends
  
  
  # Model selection
  # hessian
  if (verbose)
    cat("  Hessian\n")
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  #   err <- matrix(0, K, S)
  logDet <- 0
  
  for (m in 1:M) {
    for (k in 1:K) {
      if (k > 1) #Why start adding from index k=2???? Cause e.g. weight \pi_1 depends on the rest of the weights (they sum up to one)
        logDet <- logDet + 2.0 * log(N) - log(weights[k]) #unnormalized weights, why this added? This is the log determinant of the first block in hessian, terms corresponding to \pi_k
      #The computation of Hessian misses the regularization prior terms?
      
      hess <- hessian(lambda[[m]][k, ], Ez[k, ], binned.data[[m]], nu) #LxL matrix
      #lambda[[m]][k, ], 1xS
      #Ez[k, ], 1x1000
      #binned.data[[m]], 1000xS
      #nu 1
      ## det(H)=det(L)det(U). LU is the lower-upper decomposition of H. L is a lower triangular matrix, U is upper triangular matrix
      ##The determinant of a lower triangular matrix (or an upper triangular matrix) is the product of the diagonal entries.
      luhess <- Matrix::lu(hess) #LU decomposition
      invHess <- Matrix::solve(luhess) #inverse of Hessian matrix
      #       err[k, ] <- diag(invHess)
      #L has only ones in the diagonal so its determinant not computed
      #Why absolute, the determinant can be also negative?    
      #diagonal elements of U matrix
      #should this be det(as.matrix(Matrix::expand(luhess)$P))*Matrix::diag( Matrix::expand(luhess)$U)
      #Is the final determinant always positive, should be if we are in the extreme value
      #if the determinant is negative, this is a saddle point
      #determinant can not be zero ->
      logDet <- logDet + sum( log( abs( Matrix::diag( Matrix::expand(luhess)$U ) ) ) )
    }
  }
  
  P <- K*sum(S)+K-1 #k=S x K x M + (K-1) This should change for different every M_k?
  #gof.laplace is -log p(X|M_k) 
  gof.laplace <- last.nll + 0.5 * logDet - 0.5 * P * log(2.0 * pi); ## last.nll given by neg_log_likelihood, this is approx. -log(X|M_k)
  gof.BIC <- last.nll + 0.5 * log(N) * P #this is -BIC
  gof.AIC <- last.nll + P #this is 0.5*AIC
  gof <- c(NLE=last.nll, LogDet=logDet, Laplace=gof.laplace, BIC=gof.BIC, AIC=gof.AIC)  #goodness of fit
  
  result <- list()
  
  result$GoodnessOfFit <- gof
  result$Group <- t(Ez)
  #mixture_list <- mixture_output(binned.data, weights, lambda, err)
  #result$Mixture <- list(Weight=mixture_list$Mixture)
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
  
  #result$Fit <- list(Estimate=t(mixture_list$Estimate),
  #                   Upper=t(mixture_list$Upper),
  #                   Lower=t(mixture_list$Lower))
  
  result$Fit <- list(Estimate=lapply(lambda, function(x)t(exp(x)))) #alpha parameters
  result$Data <- binned.data #shifted and flipped data
}
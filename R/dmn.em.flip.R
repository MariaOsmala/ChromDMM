

dmn.em.flip <- function(kmeans.res, Wx, bin.width, zeta, alpha, M, K, Lx,  N, verbose, 
                   maxNumOptIter, binned.data, eta, nu, etah, nuh, numOptRelTol,  
                   EM.maxit, EM.threshold, method="BFGS", optim.options=NULL, hessian=FALSE) {
  
  
  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
  
  Ez2 <- array(0,dim=c(K,2,N))#This is 4D for each i
  
  for(i in 1:N){
    for(k in 1:K){
      print(Ez[k,i] )
      print(zeta*Ez[k,i] )
      Ez2[k,,i]=zeta*Ez[k,i]
      Ez2[k,,i]=Ez2[k,,i]/sum(Ez2[k,,i])
    }
  }
  Ez=Ez2

  Ez_list=list()
  Ez_list[[1]]=data.table::copy(Ez)
  
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
  lambda_optim_message<-vector(mode = 'list', M)
  alpha_list<-vector(mode = 'list', M)
  
  #Initialize the lambda values by optimizing Q wrt the lambda values
  for (m in 1:M) { #over data types
    lambda_optim_message[[m]]<-list()
    alpha_list[[m]]<-list()
    for (k in 1:K) { #over clusters
      if (verbose)
        setTxtProgressBar(pb, (m-1)*K+k)
      hkm <- vector(mode = 'list', maxNumOptIter+1)
      hkm_lb <- vector(mode = 'list', maxNumOptIter+1) #comes from function computation
      gradient<- vector(mode = 'list', maxNumOptIter+1)
      lb <- vector(mode = 'list', maxNumOptIter+1)
      lambda_iter <- vector(mode = 'list', maxNumOptIter+1)
      #the lambda_{kj}^{(m)} can be optimized for each k and m separately
        optim.result <- optimise_lambda_k_flip(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,,],
                                          hkm=hkm,
                                          hkm_lb=hkm_lb,
                                          gradient=gradient,
                                          lb=lb,
                                          lambda_iter=lambda_iter,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol, hessian=hessian, 
                                          method=method, optim.options=optim.options)
    
        alpha_list[[m]][[k]] =exp( do.call(rbind,lambda_iter))
        hkm <- unlist(hkm)
        hkm_lb <- unlist(hkm_lb)
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
        lambda_optim_message[[m]][[k]]<-optim.result
      
    }
  } #initialization ends
  
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  if (verbose)
    cat('\nExpectation Maximization\n')
  
  #EM loop
  iter <- 0
  last.nll <- 0
  last.nLB <- .Machine$double.xmax
  
  nll.change <- .Machine$double.xmax #1.797693e+308
  nLB.change <- .Machine$double.xmax #1.797693e+308
    
  EM_lambda_optim_message <- vector(mode = 'list', EM.maxit+1)
  EM_lambda_optim_message[[1]]<- lambda_optim_message
  
  EM_alpha_list <- vector(mode = 'list', EM.maxit+1)
  EM_alpha_list[[1]] <- alpha_list
  
  nll_list<- vector(mode = 'list', EM.maxit+1)
  nLB_list<- vector(mode = 'list', EM.maxit+1)
  
  #function neg_log_likelihood_flip is used to check the convergence
  nll_iter=1
  nll_list[[nll_iter]] <- neg_log_likelihood_flip(weights, zeta, lambda, binned.data, eta, nu, etah, nuh) #DONE!
  last.nll <- nll_list[[nll_iter]]
  nll_iter=nll_iter+1
  
  nLB_iter=1
  nLB_list[[nLB_iter]] <- neg_lower_bound_flip(Ez, weights, zeta, lambda, binned.data, eta, nu, etah, nuh)
  last.nLB <- nLB_list[[nLB_iter]]
  nLB_iter=nLB_iter+1
  
  
  nll.change_list<- vector(mode = 'list', EM.maxit+1)
  nLB.change_list<- vector(mode = 'list', EM.maxit+1)
  
  #E-step
  #M-step
  #shifting and flipping, one needs to write separate EM-loops for (shift=false, flip=false), 
  #(shift=true, flip=false), (shift=false, flip=true), (shift=true, flip=true)
  #4 different em-functions
  #while ( iter < EM.maxit ) {
   while ((iter < EM.maxit) && (nll.change > EM.threshold) && (nLB.change > EM.threshold)) {
    
    if (verbose)
      cat('Calculating Ez values...\n')
    
    #Computes the E-step, the posterior probabilities of the cluster labels
    #The old Ez values not really used
    # returns K x N matrix Z
    #why offset subtracted?
    Ez <- calc_z_flip(Ez, binned.data, weights, zeta, lambda)
    Ez_list[[iter+2]]=data.table::copy(Ez)
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
    alpha_list<-vector(mode = 'list', M)
    for (m in seq_len(M)) {
      lambda_optim_message[[m]]<-list()
      alpha_list[[m]]<-list()
      for (k in seq_len(K)) {
        if (verbose)
          setTxtProgressBar(pb, (m-1)*K+k)
        
        hkm <- vector(mode = 'list', maxNumOptIter+1)
        hkm_lb <- vector(mode = 'list', maxNumOptIter+1) #comes from function computation
        gradient <- vector(mode = 'list', maxNumOptIter+1)
        lb <- vector(mode = 'list', maxNumOptIter+1)
        lambda_iter <- vector(mode = 'list', maxNumOptIter+1)
        
        optim.result <- optimise_lambda_k_flip(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,,],
                                          hkm=hkm,
                                          hkm_lb=hkm_lb,
                                          gradient=gradient,
                                          lb=lb,
                                          lambda_iter=lambda_iter,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol, hessian=hessian, 
                                          method=method, optim.options=optim.options)
        
        alpha_list[[m]][[k]]=exp( do.call(rbind,lambda_iter))
        lambda[[m]][k,] <-optim.result$par
        lambda_optim_message[[m]][[k]]<-optim.result
        
        hkm <- unlist(hkm)
        hkm_lb <- unlist(hkm_lb)
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
    #save.image("/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/DMM-private-master-devel-works/shif.flip.debugging.RData")
    EM_lambda_optim_message[[iter+2]]<- lambda_optim_message
    EM_alpha_list[[iter+2]] <- alpha_list
    
    stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
    
    if (verbose) {
      cat('\nRead flipping state frequencies: \n')
      if( k>1 ){
        ind=apply(Ez,3,which.max)
        f = floor((ind-1) / nrow(Ez[,,1])) + 1
        print(table(f))
      }else{ k=1
      f=apply(Ez,3,which.max)
      print(table(f))
      }
      
    }
    
    weights <- rowSums(Ez) # \pi_k values Ez is K x N matrix, are these normalized? Should this be rowMeans?
    
    if (verbose)
      cat('\nCalculating negative log likelihood...\n')
    
    nll <- neg_log_likelihood_flip(weights, zeta, lambda, binned.data, eta, nu, etah, nuh)
    stopifnot(!is.na(nll), !is.infinite(nll))
    nll_list[[nll_iter]]=nll
    nll_iter=nll_iter+1
    
    nLB <- neg_lower_bound_flip(Ez, weights, zeta, lambda,
                                 binned.data, eta, nu,
                                 etah, nuh)
    
    stopifnot(!is.na(nLB), !is.infinite(nLB))
    nLB_list[[nLB_iter]]=nLB
    nLB_iter=nLB_iter+1
    
   
    nll.change <- abs(last.nll-nll)
    nLB.change <- abs(last.nLB-nLB)
    
    nll.change_list[[iter+1]]<-nll.change
    nLB.change_list[[iter+1]]<-nLB.change
    
    last.nll <- nll
    last.nLB <- nLB
    
    iter <- iter+1
    
    if (verbose)
      print(paste('--> EM Iteration:', iter, 'Neg.LL change (absolute):', round(nll.change, -log10(numOptRelTol)) ))
    print(paste('--> EM Iteration:', iter, 'Neg.LB change (absolute):', round(nLB.change, -log10(numOptRelTol)) ))
    print(paste0("Neg.LL ",round(nll,-log(numOptRelTol)) ))
    print(paste0("Neg.LL ",round(last.nLB, -log(numOptRelTol)) ))
    #  print(paste('--> EM Iteration:', iter, 'Neg.LL change:', round(nll.change, 6)))
  } #EM loop ends
  
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  nLB.data<-data.frame(iter=which(sapply(nLB_list, length)!=0), nLB=unlist(nLB_list[which(sapply(nLB_list, length)!=0)]))
  
  # Model selection
  # hessian
  #if (verbose)
  #  cat("  Hessian\n")
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  #   err <- matrix(0, K, S)
  # logDet <- 0
  # 
  # for (m in 1:M) {
  #   for (k in 1:K) {
  #     if (k > 1) #Why start adding from index k=2???? Cause e.g. weight \pi_1 depends on the rest of the weights (they sum up to one)
  #       logDet <- logDet + 2.0 * log(N) - log(weights[k]) #unnormalized weights, why this added? This is the log determinant of the first block in hessian, terms corresponding to \pi_k
  #     #The computation of Hessian misses the regularization prior terms?
  #     
  #     hess <- hessian(lambda[[m]][k, ], Ez[k, ], binned.data[[m]], nu) #LxL matrix
  #     #lambda[[m]][k, ], 1xS
  #     #Ez[k, ], 1x1000
  #     #binned.data[[m]], 1000xS
  #     #nu 1
  #     ## det(H)=det(L)det(U). LU is the lower-upper decomposition of H. L is a lower triangular matrix, U is upper triangular matrix
  #     ##The determinant of a lower triangular matrix (or an upper triangular matrix) is the product of the diagonal entries.
  #     luhess <- Matrix::lu(hess) #LU decomposition
  #     invHess <- Matrix::solve(luhess) #inverse of Hessian matrix
  #     #       err[k, ] <- diag(invHess)
  #     #L has only ones in the diagonal so its determinant not computed
  #     #Why absolute, the determinant can be also negative?    
  #     #diagonal elements of U matrix
  #     #should this be det(as.matrix(Matrix::expand(luhess)$P))*Matrix::diag( Matrix::expand(luhess)$U)
  #     #Is the final determinant always positive, should be if we are in the extreme value
  #     #if the determinant is negative, this is a saddle point
  #     #determinant can not be zero ->
  #     logDet <- logDet + sum( log( abs( Matrix::diag( Matrix::expand(luhess)$U ) ) ) )
  #   }
  # }
  
  P <- K*sum(Lx)+K-1 #k=S x K x M + (K-1) This should change for different every M_k?
  #gof.laplace is -log p(X|M_k) 
  gof.laplace <- 0 #last.nll + 0.5 * logDet - 0.5 * P * log(2.0 * pi); ## last.nll given by neg_log_likelihood, this is approx. -log(X|M_k)
  gof.BIC <- last.nll + 0.5 * log(N) * P #this is -BIC
  gof.AIC <- last.nll + P #this is 0.5*AIC
  gof <- c(NLE=last.nll, LogDet=0, Laplace=gof.laplace, BIC=gof.BIC, AIC=gof.AIC)  #goodness of fit
  
  result <- list()
  
  result$GoodnessOfFit <- gof
  result$Ez <- Ez
  result$Ez_list <- Ez_list
  result$Group <- t(apply(Ez,c(1,3),sum))
 
  result$Mixture <- list(Weight=weights/N)
  print(result$Mixture$Weight)
  result$nll.data=nll.data
  result$nLB.data=nLB.data
  result$EM_lambda_optim_message=EM_lambda_optim_message
  result$EM_alpha_list=EM_alpha_list
  
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
  
  #add flipping information
  ind=apply(Ez,3,which.max)
  if(K==1){
    s=ind
  }else{
    k = ((ind-1) %% nrow(Ez[,,1])) + 1 #this is the same as result$Group
    s = floor((ind-1) / nrow(Ez[,,1])) + 1
  }
  learned.flip.states=s
  
  unflipped.binned.data=binned.data
  for(m in 1:M){
    
    unflipped.binned.data[[m]]=matrix(0, nrow=N, ncol=Lx)
    for(i in 1:N){
      if(learned.flip.states[i]==1){
      unflipped.binned.data[[m]][i,]=binned.data[[m]][i,]
      }else{
        unflipped.binned.data[[m]][i,]=rev(binned.data[[m]][i,])
      }
    }
    rownames(unflipped.binned.data[[m]]) <- paste0('loci', seq_len(N))
    colnames(unflipped.binned.data[[m]]) <- paste0('pos', seq_len(Lx))
  }
  
  attr(unflipped.binned.data, 'flips') <- learned.flip.states
  
  
  result$Data <- unflipped.binned.data #shifted and flipped data
  
 
  result
}
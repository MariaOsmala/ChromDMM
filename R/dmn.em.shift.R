dmn.em.shift <- function(kmeans.res,  Wx, bin.width, S, xi, alpha, M, K, Lx,  N, verbose, 
                   maxNumOptIter, binned.data, eta, nu, etah, nuh, numOptRelTol, 
                   EM.maxit, EM.threshold, method="BFGS", hessian=FALSE ) {
  
  #options(digits=-log10(numOptRelTol))
  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
  #shift=true flip=false
  
  La=Wx/bin.width+S-1
  #extend the length of alpha to La
  extend=floor(S/2)
  alpha_extend<-alpha
  for(m in 1:M){
    #print(m)
    tmp=matrix(0, nrow=K, ncol=extend*2+Lx[m])
    
    # if(K==1){
    #     tmp=c(rep(0, extend), alpha[[m]], rep(0, extend))
    #     minval=min(alpha[[m]][1],alpha[[m]][ncol(alpha[[m]])])
    #     tmp[1:extend]=minval
    #     tmp[ncol(tmp)-((extend-1):0)]=minval
    #     
    #     alpha_extend[[m]]=t(as.matrix(tmp/sum(tmp))) #normalize
    #   
    # }else{ #K>1
      for(k in 1:K){
        #print(k)
        tmp[k,]=c(rep(0, extend), alpha[[m]][k,], rep(0, extend))
        minval=min(alpha[[m]][k,1],alpha[[m]][k,ncol(alpha[[m]])])
        tmp[k,1:extend]=minval
        tmp[k,ncol(tmp)-((extend-1):0)]=minval
      }
      alpha_extend[[m]]=tmp/rowSums(tmp) #normalize
    # }
    
   
    
  }
  
  alpha=alpha_extend

  Ez2 <- array(0,dim=c(K,S,N))#This is 3D for each i
  for(i in 1:N){
    for(k in 1:K){
      Ez2[k,,i]=xi*Ez[k,i]
    }
  }
  Ez=Ez2
  
  Ez_list=list()
  Ez_list[[1]]=data.table::copy(Ez)

  print(paste0("Ez_list_length: ",length(Ez_list) ))
  print(paste0("Ez_list[[1]]: ",str(Ez_list[[1]]) ))
  
  print(tracemem(Ez_list[[1]])==tracemem(Ez))
  
  
  #lambda is log(alpha)
  alpha <- lapply(alpha, function(a){a[a <= 0] <- 1e-6;a})
  lambda <- lapply(alpha, log) #natural logarithm by default
  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
  
  
  weights <- rowSums(Ez) #initial values for \bm{\pi}
  #weights <- (weights/sum(weights)*100) #these need to be numbers between 0 and zero
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
      gradient<- vector(mode = 'list', maxNumOptIter+1)
      lb <- vector(mode = 'list', maxNumOptIter+1)
      lambda_iter <- vector(mode = 'list', maxNumOptIter+1)
      #the lambda_{kj}^{(m)} can be optimized for each k and m separately
     
        optim.result=optimise_lambda_k_shift(LambdaK=lambda[[m]][k,],
                                     data=binned.data[[m]],
                                     Z=Ez[k,,],
                                     hkm=hkm,
                                     gradient=gradient,
                                     lb=lb,
                                     lambda_iter=lambda_iter,
                                     eta=eta,
                                     nu=nu,
                                     etah=etah,
                                     nuh=nuh,
                                     verbose=verbose,
                                     MAX_GRAD_ITER=maxNumOptIter,
                                     reltol=numOptRelTol, hessian=hessian, method=method)
      
        alpha_list[[m]][[k]] =exp( do.call(rbind,lambda_iter))
        hkm <- unlist(hkm)
        gradient <- unlist(gradient)
        print(length(hkm))
        print(length(gradient))
        hkm <- data.frame(Datatype=names(binned.data)[m],
                          Component=k,
                          EM.iter=0,
                          hkm=hkm, gradient=gradient, nll=optim.result$value, detH=det(optim.result$hessian))
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
  last.nll <- .Machine$double.xmax
  last.nLB <- .Machine$double.xmax
  
  nll.change <- .Machine$double.xmax #1.797693e+308
  nLB.change <- .Machine$double.xmax #1.797693e+308
  
  real.change<- .Machine$double.xmax #1.797693e+308
  nLB.real.change <- .Machine$double.xmax #1.797693e+308 
  
  
  EM_lambda_optim_message <- vector(mode = 'list', EM.maxit+1)
  EM_lambda_optim_message[[1]]<- lambda_optim_message

  EM_alpha_list <- vector(mode = 'list', EM.maxit+1)
  EM_alpha_list[[1]] <- alpha_list
  
  #use negative lower bound to check the convergence
  nLB_list<- vector(mode = 'list', EM.maxit+1)
  nll_list<- vector(mode = 'list', EM.maxit+1) #save also this, does not necessarily decrease
  
  #
  nll_iter=1
  nll_list[[nll_iter]] <- neg_log_likelihood_shift(weights, xi, lambda, binned.data, eta, nu, etah, nuh,S) #DONE!
  last.nll <- nll_list[[nll_iter]]
  nll_iter=nll_iter+1
  
  nLB_iter=1
  nLB_list[[nLB_iter]] <- neg_lower_bound_shift(Ez, weights, xi, lambda, binned.data, eta, nu, etah, nuh, S)
  last.nLB <- nLB_list[[nLB_iter]]
  nLB_iter=nLB_iter+1
  
  nll.change_list<- vector(mode = 'list', EM.maxit+1)
  nLB.change_list<- vector(mode = 'list', EM.maxit+1)
  
  #E-step
  #M-step
  #shifting and flipping, one needs to write separate EM-loops for (shift=false, flip=false), 
  #(shift=true, flip=false), (shift=false, flip=true), (shift=true, flip=true)
  #4 different em-functions
  while (iter < EM.maxit ) {
  #while ((iter < EM.maxit) && (nLB.real.change > EM.threshold) && (real.change > EM.threshold)) {
    print(paste0("EM: ",iter))
    if (verbose)
      cat('Calculating Ez values...\n')
    
    #Computes the E-step, the posterior probabilities of the cluster labels
    #The old Ez values not really used
    # returns K x N matrix Z
    #why offset subtracted?
    Ez <- calc_z_shift(Ez, binned.data, weights, xi, lambda) #TODO
    Ez_list[[iter+2]]=data.table::copy(Ez)
    
    print(paste0("Ez_list_length: ",length(Ez_list) ))
    print(paste0("E_list[[1]]: ",str(Ez_list[[1]]) ))
    
    #print(tracemem(Ez_list[[1]])==tracemem(Ez))
    
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
      #print(m)
      lambda_optim_message[[m]]<-list()
      alpha_list[[m]]<-list()
      for (k in seq_len(K)) {
        #print(k)
        if (verbose)
          setTxtProgressBar(pb, (m-1)*K+k)
        
        hkm <- vector(mode = 'list', maxNumOptIter+1)
        gradient <- vector(mode = 'list', maxNumOptIter+1)
        lb <- vector(mode = 'list', maxNumOptIter+1)
        lambda_iter <- vector(mode = 'list', maxNumOptIter+1)
        optim.result <- optimise_lambda_k_shift(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,,],
                                          hkm=hkm,
                                          gradient=gradient,
                                          lb=lb,lambda_iter=lambda_iter,
                                          eta=eta,
                                          nu=nu,
                                          etah=etah,
                                          nuh=nuh,
                                          verbose=verbose,
                                          MAX_GRAD_ITER=maxNumOptIter,
                                          reltol=numOptRelTol, hessian=hessian, method=method)
        
        alpha_list[[m]][[k]]=exp( do.call(rbind,lambda_iter))
        lambda[[m]][k,] <-optim.result$par
        lambda_optim_message[[m]][[k]]<-optim.result
        
        
        hkm <- unlist(hkm)
        
        gradient <- unlist(gradient)
        print(length(hkm))
        print(length(gradient))
        hkm <- data.frame(Datatype=names(binned.data)[m],
                          Component=k,
                          EM.iter=iter+1,
                          hkm=hkm,gradient=gradient, nll=optim.result$value, detH=det(optim.result$hessian))
        EM.diagnostics <- rbind(EM.diagnostics, hkm)
        
                
        
      }
    }
    #save.image("/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/DMM-private-master-devel-works/shif.flip.debugging.RData")
    EM_lambda_optim_message[[iter+2]]<- lambda_optim_message
    EM_alpha_list[[iter+2]] <- alpha_list

    stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))
    
   if (verbose) {
        cat('\nRead shifting state frequencies: \n')
        if( k>1 ){
            ind=apply(Ez,3,which.max)
            s = floor((ind-1) / nrow(Ez[,,1])) + 1
            print(table(s))
        }else{ k=1
          s=apply(Ez,3,which.max)
          print(table(s))
        }

    }
     
    
    
    weights <- rowSums(Ez) # \pi_k values Ez is K xSz  N matrix, are these normalized? Should this be rowMeans?
    #weights <- weights/sum(weights)
    if (verbose)
      cat('\nCalculating negative log likelihood...\n')
    
    #calculates neg. unnormalized posterior for convergence
    # weights 1xK (unnormalized, are normalized by this function/div by N)
    #lambda list of K times S matrices
    #binned.data N times S
    
   
    nll <- neg_log_likelihood_shift(weights, xi, lambda, binned.data, eta, nu, etah, nuh,S)
    stopifnot(!is.na(nll), !is.infinite(nll))
    nll_list[[nll_iter]]=nll
    nll_iter=nll_iter+1
   
    nLB <- neg_lower_bound_shift(Ez, weights, xi, lambda,
                                 binned.data, eta, nu,
                                 etah, nuh, S)
    stopifnot(!is.na(nLB), !is.infinite(nLB))
    nLB_list[[nLB_iter]]=nLB
    nLB_iter=nLB_iter+1

    
    real.change=last.nll-nll #This needs to be positive
    nLB.real.change=last.nLB-nLB
    
    if(real.change < -numOptRelTol){
      print("Warning: Neg.LL does not decrease!!!")
    }
    
    if(nLB.real.change < -numOptRelTol){
      print("Warning: Neg.LB does not decrease!!!")
    }
    
    nll.change <- abs(real.change)
    nLB.change <- abs(nLB.real.change)
    
    nll.change_list[[iter+1]]<-nll.change
    nLB.change_list[[iter+1]]<-nLB.change
    
    last.nll <- nll
    last.nLB <- nLB
    iter <- iter+1
    
    if (verbose){
      
      print(paste('--> EM Iteration:', iter, 'Neg.LL change (absolute):', round(nll.change, -log10(numOptRelTol)) ))
      print(paste('--> EM Iteration:', iter, 'Neg.LL change (real):', round(real.change, -log10(numOptRelTol)) ))
      print(paste('--> EM Iteration:', iter, 'Neg.LB change (absolute):', round(nLB.change, -log10(numOptRelTol)) ))
      print(paste('--> EM Iteration:', iter, 'Neg.LB change (real):', round(nLB.real.change, -log10(numOptRelTol)) ))
      print(paste0("Neg.LL ",round(nll,-log(numOptRelTol)) ))
      print(paste0("Neg.LL ",round(last.nLB, -log(numOptRelTol)) ))
    }
  } #EM loop ends
  print("EM finished!")
  
  # Model selection
  # hessian
  #if (verbose)
  #  cat("  Hessian\n")
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  nLB.data<-data.frame(iter=which(sapply(nLB_list, length)!=0), nLB=unlist(nLB_list[which(sapply(nLB_list, length)!=0)]))
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
  
   P <- K*sum(La)+K-1 #k=S x K x M + (K-1) This should change for different every M_k?
  #gof.laplace is -log p(X|M_k) 
  # #gof.laplace <- last.nll + 0.5 * logDet - 0.5 * P * log(2.0 * pi); ## last.nll given by neg_log_likelihood, this is approx. -log(X|M_k)
   gof.laplace <- 0
   gof.BIC <- last.nll + 0.5 * log(N) * P #this is -BIC
   gof.AIC <- last.nll + P #this is 0.5*AIC
   gof <- c(NLE=last.nll, Laplace=gof.laplace,BIC=gof.BIC, AIC=gof.AIC)  #goodness of fit remove LogDet=logDet,Laplace=gof.laplace,
  
  result <- list()
  
  result$GoodnessOfFit <- gof
  result$Ez <- Ez
  print(paste0("Ez_list_length: ",length(Ez_list) ))
  print(paste0("E_list[[1]]: ",str(Ez_list[[1]]) ))
  result$Ez_list <- Ez_list
  result$Group <- t(apply(Ez,c(1,3),sum))
  #result$Group2 <-  #sum to 1 for each n=1,..,N
  result$nll.data=nll.data
  result$nLB.data=nLB.data
  result$EM_lambda_optim_message=EM_lambda_optim_message
  result$EM_alpha_list=EM_alpha_list
  #mixture_list <- mixture_output(binned.data, weights, lambda, err)
  #result$Mixture <- list(Weight=mixture_list$Mixture)
  result$Mixture <- list(Weight=weights/N) #sum(EZ)=N
  
  EM.diagnostics <- plyr::ddply(EM.diagnostics, c('Datatype', 'Component', 'EM.iter'),
                                transform, NO.iter.count=length(hkm), NO.iter=seq_along(hkm))
  result$EM.diagnostics <- EM.diagnostics
  
  
  
  
  if(verbose) {
    print('Mixture weights: ')
    print(weights)
    print('Hard labels:')
    print(table(apply(result$Group, 1, which.max)))
    
    
    #print(table(k))
    #print(table(s))
    
  }
  
  #result$Fit <- list(Estimate=t(mixture_list$Estimate),
  #                   Upper=t(mixture_list$Upper),
  #                   Lower=t(mixture_list$Lower))
  
  result$Fit <- list(Estimate=lapply(lambda, function(x)t(exp(x)))) #alpha parameters
  
  #add shifting information
  ind=apply(Ez,3,which.max)
  if(K==1){
    s=ind
  }else{
     k = ((ind-1) %% nrow(Ez[,,1])) + 1 #this is the same as result$Group
     s = floor((ind-1) / nrow(Ez[,,1])) + 1
  }
  #result$Group_ind=k 
  #result$Shift_ind=s
  
  shift.vector=seq(-floor(S/2),floor(S/2),1)*bin.width
  learned.shift.amounts=shift.vector[s]
  
  unshifting_window=Lx[[1]]-S+1 #30
  
  unshifted.binned.data=binned.data
  
  for(m in 1:M){
    
    unshifted.binned.data[[m]]=matrix(0, nrow=N, ncol=unshifting_window)
    for(i in 1:N){
        
        unshifted.binned.data[[m]][i,]=binned.data[[m]][i, seq((S-s[i]+1) ,(Lx[[m]]-s[i]+1),1)]
    }
    rownames(unshifted.binned.data[[m]]) <- paste0('loci', seq_len(N))
    colnames(unshifted.binned.data[[m]]) <- paste0('pos', seq_len(unshifting_window))
  }
  
  attr(unshifted.binned.data, 'shifts') <- s
  
  result$Data <- unshifted.binned.data #shifted and ata
  print("dmn.em.shift finished")
  result
}

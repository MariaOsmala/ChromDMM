dmn.em.shift <- function(kmeans.res,  Wx, bin.width, S, xi, alpha, M, K, Lx,  N, verbose, 
                   maxNumOptIter, binned.data, eta, nu, etah, nuh, numOptRelTol, 
                   EM.maxit, EM.threshold, hessian=FALSE, optim.options=NULL ) {
  
  #options(digits=-log10(numOptRelTol))
  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
    
  La=Wx/bin.width+S-1
  #extend the length of alpha to La
  extend=floor(S/2)
  alpha_extend<-alpha
  for(m in 1:M){
    #print(m)
    tmp=matrix(0, nrow=K, ncol=extend*2+Lx[m])

    for(k in 1:K){
      #print(k)
      tmp[k,]=c(rep(0, extend), alpha[[m]][k,], rep(0, extend))
      minval=min(alpha[[m]][k,1],alpha[[m]][k,ncol(alpha[[m]])])
      tmp[k,1:extend]=minval
      tmp[k,ncol(tmp)-((extend-1):0)]=minval
    }
    alpha_extend[[m]]=tmp/rowSums(tmp) #normalize
  }
  
  alpha=alpha_extend

  Ez2 <- array(0,dim=c(K,S,N))#This is 3D for each i
  for(i in 1:N){
    for(k in 1:K){
      #print(dim(xi[i,]))
      Ez2[k,,i]=xi[i,]*Ez[k,i]
      Ez2[k,,i]=Ez2[k,,i]/sum(Ez2[k,,i])
    }
  }
  Ez=Ez2
  
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
     
      optim.result=optimise_lambda_k_shift(LambdaK=lambda[[m]][k,],
                                     data=binned.data[[m]],
                                     Z=Ez[k,,],
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
      
  #use negative lower bound to check the convergence
  nll_list<- vector(mode = 'list', EM.maxit+1) #save also this, does not necessarily decrease
  
  #
  nll_iter=1
  nll_list[[nll_iter]] <- neg_log_likelihood_shift(weights, xi, lambda, binned.data, eta, nu, etah, nuh,S) #DONE!
  last.nll <- nll_list[[nll_iter]]
  nll_iter=nll_iter+1
  
    
  #E-step
  #M-step
  
  while ((iter < EM.maxit) && (real.change > EM.threshold) ) {
    if (verbose)
      cat('Calculating Ez values...\n')
    
    #Computes the E-step, the posterior probabilities of the cluster labels
    Ez <- calc_z_shift(Ez, binned.data, weights, xi, lambda) #TODO
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
        optim.result <- optimise_lambda_k_shift(LambdaK=lambda[[m]][k,],
                                          data=binned.data[[m]],
                                          Z=Ez[k,,],
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
   
    real.change=last.nll-nll #This needs to be positive
    
    if(real.change < -numOptRelTol){
      print("Warning: Neg.LL does not decrease!!!")
    }
    
    last.nll <- nll
    iter <- iter+1
    
    if (verbose){
      print(paste('--> EM Iteration:', iter, 'Neg.LL change (real):', round(real.change, -log10(numOptRelTol)) ))
      print(paste0("Neg.LL ",round(nll,-log10(numOptRelTol)) ))
    }
  } #EM loop ends
    
  
  nll.data<-data.frame(iter=which(sapply(nll_list, length)!=0), nll=unlist(nll_list[which(sapply(nll_list, length)!=0)]))
  
 
  
   P <- K*sum(La)+K-1 #k=La x K x M + (K-1) This should change for different every M_k?
   gof.BIC <- last.nll + 0.5 * log(N) * P #this is -BIC
   gof.AIC <- last.nll + P #this is 0.5*AIC
   gof <- c(NLE=last.nll, BIC=gof.BIC, AIC=gof.AIC)  
  
  result <- list()
  
  result$GoodnessOfFit <- gof
  result$Ez <- Ez
  result$Group <- t(apply(Ez,c(1,3),sum))
  #result$Group2 <-  #sum to 1 for each n=1,..,N
  result$nll.data=nll.data
  result$Mixture <- list(Weight=weights/N) #sum(EZ)=N
  
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
  
  #add shifting information
  cl=apply(results$Group, 1, which.max)
  cl_ind=list()
  for(k in 1:K){
    cl_ind[[k]]=which(cl==k)
  }
  s=rep(0, length(cl))
  for(k in 1:K){
    s[ cl_ind[[k]] ]=apply(result$Ez[k,,cl_ind[[k]] ], 2, which.max ) #Ez is KxSxN
  }
  
  
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
  
  result$Data <- unshifted.binned.data #shifted and data
  result
}

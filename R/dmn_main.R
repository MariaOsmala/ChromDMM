

#' Optimize lambda without shift or flip
#' lambda_{k}^{(m)} are the lambda parameters of length S for cluster k and chromatin feature m
#'
#' @param LambdaK lambda[[m]][k,] lambda_{} current lambda values
#' @param data data=binned.data[[m]] N x S matrix
#' @param Z Ez_k, the posterior probabilities for samples originating from cluster k
#' @param hkm long list, first empty
#' @param eta
#' @param nu
#' @param etah
#' @param nuh
#' @param method default method='BFGS'
#' @param verbose
#' @param MAX_GRAD_ITER maxNumOptIter, 1000
#' @param reltol numOptRelTol 1e-12
#'
#' @return
#' @export
#'
#' @examples
#'
#'#General-purpose optimization based on Nelder–Mead,
#quasi-Newton and conjugate-gradient algorithms.
#LambdaK: initial values for the parameters to be optimized over
#fn=neg_log_evidence_lambda_pi A function to be minimized(default)/maximized,
#the first argument of the function should be the vector of parameters over which
#optimization is to take place

#gr=neg_log_derive_evidence_lambda_pi: A function to return the gradient for the BFGS
#There might be again an error with GAMMA_ITA vs (GAMMA_ITA -1),
#* GAMMA_ITA_H seems to as expected
#params,
#method=method
#control = list(maxit = MAX_GRAD_ITER, reltol = reltol)
#reltol is Relative convergence tolerance. The algorithm stops if it is unable
#to reduce the value by a factor of reltol * (abs(val) + reltol) at a step

#Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm),
#specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
#This uses function values and gradients to build up a picture of the surface to be optimized.
#a necessary condition for optimality is that the gradient be zero. not guaranteed to converge
#unless the function has a quadratic Taylor expansion near an optimum.
#However, BFGS can have acceptable performance even for non-smooth optimization instances.
#the Hessian matrix is approximated using updates specified by gradient evaluations (or approximate gradient evaluations)

#Nocedal and Wright (1999) is a comprehensive reference for the previous three methods.
#Nocedal, J. and Wright, S. J. (1999). Numerical Optimization. Springer.

#neg_log_evidence_lambda_pi:  computes the value of the negative expected log posterior
#i.e. lower bound to be optimized
# i.e. the whole expected log posterior Q. Does not consider \pi.
#Computes the value for one cluster k and datatype m. There might be errors(regularization term incorrect,
#not all alpha terms included), check
#'
optimise_lambda_k <- function(LambdaK, data, Z, hkm, gradient, eta, nu,
                              etah, nuh, method='BFGS',
                              verbose=FALSE, MAX_GRAD_ITER=1000,
                              reltol = 1e-12, optim.options=NULL, hessian=FALSE) {

  if(length(optim.options)==0){
    optim.options=c(maxit=MAX_GRAD_ITER, reltol=reltol)
  }
  
  hkm_index <- vector(mode='integer', 1)
  gradient_index <- vector(mode='integer', 1)
  params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                 hkm=hkm, gradient=gradient, hkm_index=hkm_index, gradient_index=gradient_index)
  
  fn=neg_log_evidence_lambda_pi
  gr=neg_log_derive_evidence_lambda_pi
  
  optim.result <- optim(LambdaK, fn=fn,
                        gr=gr, params,
                        method=method, control = optim.options, hessian=hessian)

  if(optim.result$convergence != 0){
    warning('!!!!! Numerical Optimization did not converge !!!!!!!\n')
    if(optim.result$convergence == 1){
      warning('!!!!! iteration limit maxit had been reached !!!!!!!\n')
    }
    
  }
  
  return(optim.result)
}

#' Optimize lambda with shift without flip
#' lambda_{k}^{(m)} are the lambda parameters of length La for cluster k and chromatin feature m
#'
#' @param LambdaK lambda[[m]][k,] lambda_{} current lambda values
#' @param data data=binned.data[[m]] N x Lx matrix
#' @param Z Ez_k, the posterior probabilities for samples originating from cluster k, Z<-array(0,dim=c(S,N))
#' @param hkm long list, first empty
#' @param eta
#' @param nu
#' @param etah
#' @param nuh
#' @param method default method='BFGS'
#' @param verbose
#' @param MAX_GRAD_ITER maxNumOptIter, 1000
#' @param reltol numOptRelTol 1e-12
#'
#' @return
#' @export
#'
#' @examples

optimise_lambda_k_shift <- function(LambdaK, data, Z, hkm, gradient, eta, nu,
                                         etah, nuh, method='BFGS',
                                         verbose=FALSE, MAX_GRAD_ITER=1000,
                                         reltol = 1e-12, hessian=FALSE, optim.options=NULL) {
  
  if(length(optim.options)==0){
    optim.options=c(maxit=MAX_GRAD_ITER, reltol=reltol)
  }
  
  hkm_index <- vector(mode='integer', 1)
  gradient_index <- vector(mode='integer', 1)
  
  params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                 hkm=hkm, gradient=gradient, hkm_index=hkm_index,
                 gradient_index=gradient_index)
  
  fn=neg_log_evidence_lambda_pi_shift
  gr=neg_log_derive_evidence_lambda_pi_shift
  
  optim.result <- optim(LambdaK,fn=fn, gr=gr, lambda_iter, params,
                        method=method, control = optim.options, hessian=hessian)

  
  if(optim.result$convergence != 0){
    warning('!!!!! Numerical Optimization did not converge !!!!!!!\n')
    if(optim.result$convergence == 1){
      warning('!!!!! iteration limit maxit had been reached !!!!!!!\n')
    }
    
  }
  
  return(optim.result)
}

#' Optimize lambda with flip without shift
#' lambda_{k}^{(m)} are the lambda parameters of length La for cluster k and chromatin feature m
#'
#' @param LambdaK lambda[[m]][k,] lambda_{} current lambda values
#' @param data data=binned.data[[m]] N x Lx matrix
#' @param Z Ez_k, the posterior probabilities for samples originating from cluster k,
#' @param hkm long list, first empty
#' @param eta
#' @param nu
#' @param etah
#' @param nuh
#' @param method default method='BFGS'
#' @param verbose
#' @param MAX_GRAD_ITER maxNumOptIter, 1000
#' @param reltol numOptRelTol 1e-12
#'
#' @return
#' @export
#'
#' @examples



optimise_lambda_k_flip <- function(LambdaK, data, Z, hkm, gradient, eta, nu,
                                    etah, nuh, method='BFGS',
                                    verbose=FALSE, MAX_GRAD_ITER=1000,
                                    reltol = 1e-12, hessian=FALSE, optim.options=NULL) {
  

  
  if(length(optim.options)==0){
    optim.options=c(maxit=MAX_GRAD_ITER, reltol=reltol)
  }
  
  hkm_index <- vector(mode='integer', 1)
  gradient_index <- vector(mode='integer', 1)
  
  params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                 hkm=hkm,  gradient=gradient, hkm_index=hkm_index,
                 gradient_index=gradient_index)
  
  fn=neg_log_evidence_lambda_pi_flip
  gr=neg_log_derive_evidence_lambda_pi_flip
  
  
  optim.result <- optim(LambdaK,fn=fn, gr=gr, lambda_iter, params,
                        method=method, control = optim.options, hessian=hessian)
  
  
  if(optim.result$convergence != 0){
    warning('!!!!! Numerical Optimization did not converge !!!!!!!\n')
    if(optim.result$convergence == 1){
      warning('!!!!! iteration limit maxit had been reached !!!!!!!\n')
    }
  }
  return(optim.result)
}


#' Optimize lambda with shift and flip
#' lambda_{k}^{(m)} are the lambda parameters of length La for cluster k and chromatin feature m
#'
#' @param LambdaK lambda[[m]][k,] lambda_{} current lambda values
#' @param data data=binned.data[[m]] N x Lx matrix
#' @param Z Ez_k, the posterior probabilities for samples originating from cluster k, Z<-array(0,dim=c(S,2,N))
#' @param hkm long list, first empty
#' @param eta
#' @param nu
#' @param etah
#' @param nuh
#' @param method default method='BFGS'
#' @param verbose
#' @param MAX_GRAD_ITER maxNumOptIter, 1000
#' @param reltol numOptRelTol 1e-12
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
optimise_lambda_k_shift_flip <- function(LambdaK, data, Z, hkm, gradient, eta, nu,
                              etah, nuh, method='BFGS',
                              verbose=FALSE, MAX_GRAD_ITER=1000,
                              reltol = 1e-12, hessian=FALSE, optim.options=NULL) {
  
  if(length(optim.options)==0){
    optim.options=c(maxit=MAX_GRAD_ITER, reltol=reltol)
  }

  hkm_index <- vector(mode='integer', 1)
  gradient_index <- vector(mode='integer', 1)
   
  params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                 hkm=hkm, gradient=gradient, hkm_index=hkm_index,
                 gradient_index=gradient_index)
  
  fn=neg_log_evidence_lambda_pi_shift_flip
  gr=neg_log_derive_evidence_lambda_pi_shift_flip
  
  optim.result <- optim(LambdaK, fn=fn,
                        gr=gr, lambda_iter, params,
                        method=method,  control = optim.options, hessian=hessian)
  if(optim.result$convergence != 0){
    warning('!!!!! Numerical Optimization did not converge !!!!!!!\n')
    if(optim.result$convergence == 1){
      warning('!!!!! iteration limit maxit had been reached !!!!!!!\n')
    }
  }
  
  return(optim.result)
}



#' DMN.cluster fit the DirichletMultinomial mixture model
#' The user of the package can not use this function directly but through dmn
#'
#' @param count signals.subset, a list of chromatin feature signals,
#' each element is a N x window matrix. If matrix, converted to a list
#' @param K the number of clusters, scalar not vector
#' @param bin.width bin_size 40 (default 50)
#' @param S the number of shift states, if shift is true, this needs to be odd number of at least 3, default 1(no shift)
#' @param verbose Print progress as "DMN, K=%d, M=%d, Iteration=%d", k, M, r) default FALSE
#' @param seed default false
#' @param shift default false
#' @param flip default false
#' @param eta default 0.1
#' @param nu default 0.1
#' @param etah default 0
#' @param nuh default 0
#' @param maxIt default 250
#' @param EM.threshold default 1e-6
#' @param soft.kmeans.maxit default 1000
#' @param soft.kmeans.stiffness default 50
#' @param randomInit default true
#' @param repetition default 4
#' @param maxNumOptIter default 1000
#' @param numOptRelTol default 1e-12
#' @param parallel default true
#' @param method default BFGS, option gradient.descent
#' @return
#' @export
#'
#' @examples


DMN.cluster <- function(count.data,
                        K,
                        bin.width=50,
                        S=1,
                        seed=F,
                        shift.reads=F,
                        flip=F,
                        parallel.shift=F,
                        verbose=F,
                        eta=0.1, nu=0.1,
                        etah=0, nuh=0,
                        xi=NULL, zeta=NULL,
                        EM.maxit=250, EM.threshold=1e-6,
                        soft.kmeans.maxit=1000, soft.kmeans.stiffness=50,
                        randomInit=T,
                        maxNumOptIter=1000, numOptRelTol=1e-12,init="random", method="BFGS", 
                        optim.options=NULL, hessian=FALSE,..., learning.rate=1e-3) {

  if (seed != F) set.seed(seed)

  if(!is.list(count.data))
    stop('count.data must be a list!')

  zero.indices <- which(apply(sapply(count.data, rowSums), 1, prod)==0)
  if (length(zero.indices) > 0) {
    warning(sprintf('Excluding are all-zero row(s) in the data: %s\n',
                    paste(names(zero.indices), collapse = ',')))


    count.data <- lapply(count.data, function(cd)cd[-zero.indices,])
  }

  if(!shift.reads) S <- 1

  if (shift.reads && ( S%%2 == 0 || S==1 ))
    stop('Number of shift states must be odd and at least 3')

  M <- length(count.data)

  if (M < 1)
    stop('List must have at least 1 element')

  if(is.null(names(count.data)))
    stop('List elements must have a name')

  if (!all(sapply(count.data, is.matrix)))
    stop('All list items must be a matrix')

  if(!all(diff(sapply(count.data, nrow)) == 0))
    stop('Row numbers of matrices are not same!')

  #disable GSL error handler, it aborts the program in case of an error
  disable_gsl_error_handler()

  N <- nrow(count.data[[1]]) #number of samples
  Wx <- sapply(count.data, ncol) #original window W_x

  if (sum(Wx %% bin.width) != 0)
    stop('Data column length is not a multiple of bin width')

  # These windows represent the portion of the data that we are
  # actually working on. Shifting the data simply means shifting these windows.
  # For example if S=1000 and shift.ratio=1/10, then inner.window
  # range is 50-950 meaning that we can shift the data in -50,+50 range
  #
  # TODO: Add support for matrices with different num. of columns
  # 50*(1/16) =3
  #left.limit <- right.limit <- floor( (min(Wx)/bin.width) * (shift.ratio/2) )
  left.limit <- right.limit <- floor( (min(Wx)/bin.width) * (0/2) )
  inner.windows <- IRangesList(start=as.list(rep(1, N)),
                               end=as.list(rep(min(Wx), N)))
  inner.windows <- resize(inner.windows,
                          min(Wx)-((left.limit + right.limit)*bin.width), #2000 - (3+3)*40 =2000-240=1760
                          fix='center')

  #initialize binned data and do row/col naming
  binned.data <- Map(function(data, ind){
    m <- matrix(0, N, width(inner.windows[[ind]])/bin.width) #N x (Wx/B=50)
    rownames(m) <- paste0('loc', seq_len(N))
    colnames(m) <- paste0('bin', seq_len(ncol(m)))
    m
  }, count.data, seq_len(M))

  #maxlimits <- left.limit * bin.width #120
  #minlimits <- -right.limit * bin.width #-120

  # add nonbinned data, windows, shifts and bin width as attributes
  # to the binned data
  attr(binned.data, 'nonbinned') <- count.data #N x window matrices
  attr(binned.data, 'windows') <- inner.windows
  attr(binned.data, 'shifts') <- numeric(N) #optimal shift state
  attr(binned.data, 'flips') <- logical(N) #boolean vector, optimal flip state
  attr(binned.data, 'bin.width') <- bin.width
  #attr(binned.data, 'shift.limits') <- matrix(c(minlimits, maxlimits),
  #                                            nrow=N, ncol=2, byrow=T)

  #bin data and define shifting function
  #extract_binned_signal converts numerical matrix into an integer matrix
  #Does nothing?
  binned.data <-extract_binned_signal(binned.data, seq_len(N)) #binned.data$H3K4me1: NxS matrix
  Lx <- sapply(binned.data, ncol)

  #row-wise normalization of all datatypes for soft-kmeans

  #mat N x 44
  #for each row, compute the sum
  #divide the values in each row by the sum
  if(init=="kmeans++"){
    kmeans.binned.data <- mapply(function(mat, name){
      t( apply(mat, 1, function(row) {
        s<-sum(row);
        if(s!=0) row/s else row}) )
      }, binned.data, names(binned.data), SIMPLIFY=F)
  
    
    #concatenate data for soft-kmeans
    kmeans.binned.data <- do.call(cbind, kmeans.binned.data)
    #K=2
    kmeanspp.centers <- kmeanspp_initialize(as.matrix(kmeans.binned.data), K) #indices of the centers
    kmeanspp.centers <- kmeans.binned.data[kmeanspp.centers, , drop=F] #the actual centers, K x (M*S)
    
    #rowNorm=F, the rows were already normalized
    kmeans.res <- soft_kmeans(kmeans.binned.data, K, verbose=verbose, #This is in dmn.cpp file
                              randomInit=randomInit, centers=kmeanspp.centers,
                              stiffness=soft.kmeans.stiffness, rowNorm=F)
    #list of 3
    #$centers
    #$weights
    #$labels K x N matrix
    rm(kmeans.binned.data)
  
    if(verbose) {
     cat('k-means hard label frequencies:')
     print(table(apply(kmeans.res$labels, 2, which.max)))
    }
    alpha <- kmeans.res$centers #K x (M*S)
    stopifnot(!is.na(alpha), !is.infinite(alpha))
    
    #split centers given by soft kmeans, i.e. unconcatenate
    col.ends <- cumsum(Lx)
    col.starts <- col.ends - Lx + 1
    alpha <- mapply(function(s,e)alpha[,s:e,drop=F], col.starts, col.ends, SIMPLIFY = F) #List of M
  
    #kmeans++ init
  }
  
  
  # if(init=="squeeze"){
  #   kmeans.binned.data <- mapply(function(mat, name){
  #     t( apply(mat, 1, function(row) {
  #       s<-sum(row);
  #       if(s!=0) row/s else row}) )
  #   }, binned.data, names(binned.data), SIMPLIFY=F)
  #   
  #   
  #   #concatenate data for soft-kmeans
  #   kmeans.binned.data <- do.call(cbind, kmeans.binned.data)
  #   #K=2
  #   kmeanspp.centers <- kmeanspp_initialize(as.matrix(kmeans.binned.data), K) #indices of the centers
  #   kmeanspp.centers <- kmeans.binned.data[kmeanspp.centers, , drop=F] #the actual centers, K x (M*S)
  #   
  #   #rowNorm=F, the rows were already normalized
  #   kmeans.res <- soft_kmeans(kmeans.binned.data, K, verbose=verbose, #This is in dmn.cpp file
  #                             randomInit=randomInit, centers=kmeanspp.centers,
  #                             stiffness=soft.kmeans.stiffness, rowNorm=F)
  #   #list of 3
  #   #$centers
  #   #$weights
  #   #$labels K x N matrix
  #   rm(kmeans.binned.data)
  #   
  #   if(verbose) {
  #     cat('k-means hard label frequencies:')
  #     print(table(apply(kmeans.res$labels, 2, which.max)))
  #   }
  #   alpha <- kmeans.res$centers #K x (M*S)
  #   stopifnot(!is.na(alpha), !is.infinite(alpha))
  #   
  #   #split centers given by soft kmeans, i.e. unconcatenate
  #   col.ends <- cumsum(Lx)
  #   col.starts <- col.ends - Lx + 1
  #   alpha <- mapply(function(s,e)alpha[,s:e,drop=F], col.starts, col.ends, SIMPLIFY = F) #List of M
  #   
  #   squeezed_alpha <- list()
  #   for(m in 1:M){
  #     squeezed_alpha[[ names(binned.data)[m] ]]=matrix(0, nrow=K, ncol=Lx[m])
  #     
  #     center.left=Lx[m]/2
  #     center.right=center.left+1
  #     for(k in 1:K){
  #       squeezed_alpha[[ names(binned.data)[m] ]][k,]=alpha[[names(binned.data)[m] ]][k,]
  #       
  #       #floor(S/2) first non-zero elements
  #       first.nonzero=which(alpha[[names(binned.data)[m]] ][k,]!=0)[1:floor(S/2)]
  #       last.nonzero=tail(which(alpha[[names(binned.data)[m]] ][k,]!=0), floor(S/2))
  #       
  #       inner.index=seq( (first.nonzero[1]+1), ( tail(last.nonzero,1)-1), 1 )
  #       
  #       value.index=as.vector(c(rep(first.nonzero[1], floor(S/2) ),
  #                               seq(first.nonzero[1]+1, center.left-1,2 ),
  #                               seq(center.right+1, tail(last.nonzero,1)-1,2 ),
  #                               rep(tail(last.nonzero,1), floor(S/2) )))
  #       
  #       squeezed_alpha[[ names(binned.data)[m] ]][k, inner.index]=alpha[[names(binned.data)[m]] ][k,value.index]
  #       squeezed_alpha[[ names(binned.data)[m] ]][k,]=squeezed_alpha[[ names(binned.data)[m] ]][k,]/sum(squeezed_alpha[[ names(binned.data)[m] ]][k,])
  #     }
  #     colnames(squeezed_alpha[[ names(binned.data)[m] ]]) <- paste0('bin', seq_len(Lx[m]))
  #   }
  #   
  #   alpha=squeezed_alpha 
  #   #squeeze init
  #   
  # }
  # 
  if(init=="uniform"){
    
    alpha <- list()
    for(m in 1:M){
      alpha[[ names(binned.data)[m] ]]=matrix(0, nrow=K, ncol=Lx[m])
      for(k in 1:K){
        alpha[[ names(binned.data)[m] ]][k,]=rep(1, Lx[m])/Lx[m]
      }
      colnames(alpha[[ names(binned.data)[m] ]]) <- paste0('bin', seq_len(Lx[m]))
    }
    
    
    
    kmeans.res <- list()
    kmeans.res$labels=matrix(0, nrow=K, ncol=N)
    for(k in 1:(K-1)){
      kmeans.res$labels[k,]=sample( seq(0.1,0.9,0.1), N, replace=TRUE)
    }
    kmeans.res$labels[K,]=1-colSums(kmeans.res$labels)
    colnames(kmeans.res$labels) <- paste0('loc', seq_len(N))
  }
  
  if(init=="random"){
    means.ind=sample.int(N, size=K, replace=TRUE)
    alpha <- list()
    for(m in 1:M){
      alpha[[ names(binned.data)[m] ]]=matrix(0, nrow=K, ncol=Lx[m])
      for(k in 1:K){
        alpha[[ names(binned.data)[m] ]][k,]=binned.data[[m]][means.ind[k],]/sum(binned.data[[m]][means.ind[k],])
      }
      colnames(alpha[[ names(binned.data)[m] ]]) <- paste0('bin', seq_len(Lx[m]))
    }
    
    kmeans.res <- list()
    kmeans.res$labels=matrix(0, nrow=K, ncol=N)
    for(k in 1:(K-1)){
      kmeans.res$labels[k,]=sample( seq(0.1,0.9,0.1), N, replace=TRUE)
    }
    kmeans.res$labels[K,]=1-colSums(kmeans.res$labels)
    colnames(kmeans.res$labels) <- paste0('loc', seq_len(N))
  }
  
  if(shift.reads==TRUE && flip==TRUE){
    print("shifting and flipping")
    result=dmn.em.shift.flip(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, S=S, xi=xi, zeta=zeta, alpha=alpha, 
                        M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                        maxNumOptIter=maxNumOptIter, binned.data=binned.data, 
                        eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                        EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method, optim.options=optim.options, hessian=hessian)
    
  }else if(shift.reads==TRUE && flip==FALSE){ #only shift
  
    if(method=="gradient.descent"){
      result=dmn.em.shift.gradient.descent(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, S=S, xi=xi, alpha=alpha, 
                          M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                          maxNumOptIter=maxNumOptIter, binned.data=binned.data, 
                          eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                          EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method,  learning.rate=learning.rate)
    }else{
      result=dmn.em.shift(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, S=S, xi=xi, alpha=alpha, 
                          M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                          maxNumOptIter=maxNumOptIter, binned.data=binned.data, 
                          eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                          EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method, optim.options=optim.options, hessian=hessian)
    }
  
    
  }else if(shift.reads==FALSE && flip==TRUE){ #only flip
    
    result=dmn.em.flip(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, zeta=zeta, alpha=alpha, 
                        M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                        maxNumOptIter=maxNumOptIter, binned.data=binned.data, 
                        eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                        EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method, optim.options=optim.options, hessian=hessian)
    
    
  }else{ # shift.reads==FALSE && flip==FALSE
    
    if(method=="gradient.descent"){
      result = dmn.em.gradient.descent(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, alpha=alpha, M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                      maxNumOptIter=maxNumOptIter, binned.data=binned.data, eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                      EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method,  learning.rate=learning.rate)
    }else{
      result = dmn.em(kmeans.res=kmeans.res, Wx=Wx, bin.width=bin.width, alpha=alpha, M=M, K=K, Lx=Lx, N=N, verbose=verbose, 
                      maxNumOptIter=maxNumOptIter, binned.data=binned.data, eta=eta, nu=nu, etah=etah, nuh=nuh, numOptRelTol=numOptRelTol, 
                      EM.maxit=EM.maxit, EM.threshold=EM.threshold, method=method, optim.options=optim.options, hessian=hessian)
    }
    
  }

  
  return(result)
}

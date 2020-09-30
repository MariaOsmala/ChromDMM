

#' Optimize lambda
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
#'
#'
optimise_lambda_k <- function(LambdaK, data, Z, hkm, eta, nu,
                              etah, nuh, method='BFGS',
                              verbose=FALSE, MAX_GRAD_ITER=1000,
                              reltol = 1e-12) {

  hkm_index <- vector(mode='integer', 1)
  params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                 hkm=hkm, hkm_index=hkm_index)
  #General-purpose optimization based on Nelder–Mead,
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
  optim.result <- optim(LambdaK, fn=neg_log_evidence_lambda_pi,
                        gr=neg_log_derive_evidence_lambda_pi, params,
                        method=method, control = list(maxit = MAX_GRAD_ITER,
                                                      reltol = reltol))
  #returns a list with components
  #par:The best set of parameters found.
  #value:The value of fn corresponding to par.

  if(optim.result$convergence != 0)
    warning('!!!!! Numerical Optimization did not converge !!!!!!!\n')

  return(optim.result$par)
}

#' DMN.cluster fit the DirichletMultinomial mixture model
#' The user of the package can not use this function directly but through dmn
#'
#' @param count signals.subset, a list of chromatin feature signals,
#' each element is a N x window matrix. If matrix, converted to a list
#' @param K the number of clusters, scalar not vector
#' @param bin.width bin_size 40 (default 50)
#' @param shift.ratio default 1/8
#' @param verbose Print progress as "DMN, K=%d, M=%d, Iteration=%d", k, M, r) default FALSE
#' @param seed default false
#' @param shift.reads default true
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
#'
#' @return
#' @export
#'
#' @examples


DMN.cluster <- function(count.data,
                        K,
                        bin.width=50,
                        shift.ratio=1/8,
                        seed=F,
                        shift.reads=T,
                        flip=F,
                        parallel.shift=F,
                        verbose=F,
                        eta=0.1, nu=0.1,
                        etah=0, nuh=0,
                        EM.maxit=250, EM.threshold=1e-6,
                        soft.kmeans.maxit=1000, soft.kmeans.stiffness=50,
                        randomInit=T,
                        maxNumOptIter=1000, numOptRelTol=1e-12) {

  if (seed != F) set.seed(seed)

  if(!is.list(count.data))
    stop('count.data must be a list!')

  zero.indices <- which(apply(sapply(count.data, rowSums), 1, prod)==0)
  if (length(zero.indices) > 0) {
    warning(sprintf('Excluding are all-zero row(s) in the data: %s\n',
                    paste(names(zero.indices), collapse = ',')))


    count.data <- lapply(count.data, function(cd)cd[-zero.indices,])
  }

  if(!shift.reads) shift.ratio <- 0

  if (shift.reads && (shift.ratio > 1 || shift.ratio < 0 ))
    stop('Shift ratio must be between 0 and 1')

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
  S.nonbinned <- sapply(count.data, ncol) #original window

  if (sum(S.nonbinned %% bin.width) != 0)
    stop('Data column length is not a multiple of bin width')

  # These windows represent the portion of the data that we are
  # actually working on. Shifting the data simply means shifting these windows.
  # For example if S=1000 and shift.ratio=1/10, then inner.window
  # range is 50-950 meaning that we can shift the data in -50,+50 range
  #
  # TODO: Add support for matrices with different num. of columns
  # 50*(1/16) =3
  left.limit <- right.limit <- floor( (min(S.nonbinned)/bin.width) * (shift.ratio/2) )

  inner.windows <- IRangesList(start=as.list(rep(1, N)),
                               end=as.list(rep(min(S.nonbinned), N)))
  inner.windows <- resize(inner.windows,
                          min(S.nonbinned)-((left.limit + right.limit)*bin.width), #2000 - (3+3)*40 =2000-240=1760
                          fix='center')

  #initialize binned data and do row/col naming
  binned.data <- Map(function(data, ind){
    m <- matrix(0, N, width(inner.windows[[ind]])/bin.width) #N x (1760/40 =44)
    rownames(m) <- paste0('loc', seq_len(N))
    colnames(m) <- paste0('bin', seq_len(ncol(m)))
    m
  }, count.data, seq_len(M))

  maxlimits <- left.limit * bin.width #120
  minlimits <- -right.limit * bin.width #-120

  # add nonbinned data, windows, shifts and bin width as attributes
  # to the binned data
  attr(binned.data, 'nonbinned') <- count.data #N x window matrices
  attr(binned.data, 'windows') <- inner.windows
  attr(binned.data, 'shifts') <- numeric(N)
  attr(binned.data, 'flips') <- logical(N) #boolean vector
  attr(binned.data, 'bin.width') <- bin.width
  attr(binned.data, 'shift.limits') <- matrix(c(minlimits, maxlimits),
                                              nrow=N, ncol=2, byrow=T)

  #bin data and define shifting function
  #extract_binned_signal converts numerical matrix into an integer matrix
  #Does nothing?
  binned.data <-extract_binned_signal(binned.data, seq_len(N)) #binned.data$H3K4me1: NxS matrix
  S <- sapply(binned.data, ncol)

  #row-wise normalization of all datatypes for soft-kmeans

  #mat N x 44
  #for each row, compute the sum
  #divide the values in each row by the sum

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
  col.ends <- cumsum(S)
  col.starts <- col.ends - S + 1
  alpha <- mapply(function(s,e)alpha[,s:e,drop=F], col.starts, col.ends, SIMPLIFY = F) #List of M

  #lambda is log(alpha)
  alpha <- lapply(alpha, function(a){a[a <= 0] <- 1e-6;a})
  lambda <- lapply(alpha, log) #natural logarithm by default

  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))

  Ez <- kmeans.res$labels #K x N initial values of the posterior probabilities of cluster assignments
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

  #Initialize the lambda values by optimizing Q wrt the lambda values
  for (m in 1:M) { #over data types
    for (k in 1:K) { #over clusters
      if (verbose)
        setTxtProgressBar(pb, (m-1)*K+k)
      hkm <- vector(mode = 'list', maxNumOptIter+1)

      #the lambda_{kj}^{(m)} can be optimized for each k and m separately
      lambda[[m]][k,] <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
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
  }

  stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))

  if (verbose)
    cat('\nExpectation Maximization\n')

  #EM loop
  iter <- 0
  last.nll <- 0
  nll.change <- .Machine$double.xmax #1.797693e+308
  EM.diagnostics <- data.frame()
  #E-step
  #M-step
  #shifting and flipping
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

    for (m in seq_len(M)) {
      for (k in seq_len(K)) {
        if (verbose)
          setTxtProgressBar(pb, (m-1)*K+k)

        hkm <- vector(mode = 'list', maxNumOptIter+1)
        lambda[[m]][k,] <- optimise_lambda_k(LambdaK=lambda[[m]][k,],
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




        hkm <- unlist(hkm)
        hkm <- data.frame(Datatype=names(binned.data)[m],
                          Component=k,
                          EM.iter=iter,
                          hkm=hkm)
        EM.diagnostics <- rbind(EM.diagnostics, hkm)
      }
    }
    #save.image("/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/DMM-private-master-devel-works/shif.flip.debugging.RData")

    stopifnot(!is.na(unlist(lambda)), !is.infinite(unlist(lambda)))

    if ((shift.ratio > 0 && shift.reads) || flip){
      if (verbose)
        cat('\nOptimizing shift/flip parameters...\n')
      #this is in util.R
      #optimization_func is in dmn.cpp, what does it try to optimize?
      #returns $shift (length N vector)
      #returns $flip (length N locigal)
      best.shift.and.flip.params <- brute.force.integer.minimizer(fn=optimization_func,
                                                     limits=attr(binned.data, 'shift.limits'),
                                                     step=bin.width,
                                                     parallel.shift=parallel.shift,
                                                     verbose=verbose,
                                                     shift.reads=shift.reads,
                                                     flip=flip,
                                                     data=binned.data,
                                                     alpha=lapply(lambda, exp),
                                                     Z=Ez)
      if (verbose) {
        cat('\nRead shifting frequencies: \n')
        print(table(best.shift.and.flip.params$shift))

        cat('\nRead flip parameters: \n')
        print(table(best.shift.and.flip.params$flip))

      }
      binned.data <- shift.and.flip.signal(binned.data, seq_len(N),
                                           best.shift.and.flip.params$shift,
                                           best.shift.and.flip.params$flip)
    } #end of shifting and flipping

    weights <- rowSums(Ez) # \pi_k values Ez is K x N matrix, are these normalized? Should this be rowMeans?

    if (verbose)
      cat('\nCalculating negative log likelihood...\n')

    #calculates neg. unnormalized posterior for convergence
    # weights 1xK (unnormalized, are normalized by this function/div by N)
    #lambda list of K times S matrices
    #binned.data N times S
    nll <- neg_log_likelihood(weights, lambda, binned.data, eta, nu, etah, nuh)
    stopifnot(!is.na(nll), !is.infinite(nll))

    nll.change <- abs(last.nll - nll)
    last.nll <- nll

    iter <- iter+1

    if (verbose)
      print(paste('--> EM Iteration:', iter, 'Neg.LL change:', round(nll.change, 6)))
  } #EM loop ends


  # Model selection
  # hessian
  if (verbose)
    cat("  Hessian\n")

#   err <- matrix(0, K, S)
  logDet <- 0

  for (m in 1:M) {
    for (k in 1:K) {
      if (k > 1)
        logDet <- logDet + 2.0 * log(N) - log(weights[k]) #unnormalized weights, should one normalize first?
      #The computation of Hessian misses the regularization prior terms?
      #Why Ez[k,] used in the computation of Hessian? What is the function to be derived twice? Energy function?
      hess <- hessian(lambda[[m]][k, ], Ez[k, ], binned.data[[m]], nu) #SxS matrix
      #lambda[[m]][k, ], 1xS
      #Ez[k, ], 1x1000
      #binned.data[[m]], 1000xS
      #nu 1
      luhess <- Matrix::lu(hess) #LU decomposition
      invHess <- Matrix::solve(luhess) #inverse of Hessian matrix
#       err[k, ] <- diag(invHess)
      logDet <- logDet + sum(log(abs(Matrix::diag(Matrix::expand(luhess)$U))))
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

  return(result)
}

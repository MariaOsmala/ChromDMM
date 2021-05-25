setClass("DMN",
    representation=representation(goodnessOfFit="numeric",
                                  group="matrix",
                                  mixture="list",
                                  fit="list",
                                  Ez="array",
                                  EM.diagnostics="data.frame",
                                  Data="list",
                                  nll.data="data.frame",
                                  nLB.data="data.frame",
                                  EM_lambda_optim_message="list",
                                  EM_alpha_list="list",
                                  Ez_list="list"))

.DMN <-
    function(goodnessOfFit, group, mixture, fit, ...)
{
    new("DMN", goodnessOfFit=goodnessOfFit,
        group=group,
        mixture=mixture,
        fit=fit,
        ...)
}

#' The main function to fit the DirichletMultinomial mixture model
#'
#' @param count signals.subset, a list of chromatin feature signals, 
#' each element is a N x window matrix. If matrix, converted to a list
#' @param K vector of the number of clusters
#' @param bin.width bin_size 40 (default 50)
#' @param S the number of shift states, if shift true, this needs to be odd number of at least 3, default 1(no shift)
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
#'@param init default "random"
#'
#' @return
#' @export
#'
#' @examples


dmn <-
    function(count,
             K,
             bin.width=50,
             S=1,
             verbose=FALSE,
             seed=F,
             shift.reads=F,
             flip=F,
             eta=0.1,
             nu=0.1,
             etah=0,
             nuh=0,
             xi=NULL,
             maxIt=250,
             EM.threshold=1e-6,
             soft.kmeans.maxit=1000,
             soft.kmeans.stiffness=50,
             randomInit=T,
             repetition=4,
             maxNumOptIter=1000,
             numOptRelTol=1e-12,
             parallel=T, init="random", method="BFGS", hessian=FALSE, learning.rate=1e-3)
{
    if (is.matrix(count)) count <- list(Data=count)
    stopifnot(is.list(count))
    
    M <- length(count) #number of chromatin features

    if ( (randomInit==F) || seed || ( (length(K) == 1) && (K==1) ) ){
      repetition<-1
    }
    repet.func <- function(r, k, ...){
      if (verbose){
        message(sprintf("DMN, K=%d, M=%d, Iteration=%d", k, M, r))
      }

      for (i in seq(M)) {
        if (mode(count[[i]]) == 'numeric') { 
          mode(count[[i]]) <- "integer" #set the type or storage mode of an object
        }
      }
      
      ans <- DMN.cluster(count.data=count,
                         K=k,
                         bin.width=bin.width,
                         S=S,
                         seed=seed,
                         shift.reads=shift.reads,
                         flip=flip,
                         verbose=verbose,
                         eta=eta, nu=nu,
                         etah=etah, nuh=nuh, xi=xi,
                         EM.maxit=maxIt, EM.threshold=EM.threshold,
                         soft.kmeans.maxit=soft.kmeans.maxit,
                         soft.kmeans.stiffness=soft.kmeans.stiffness,
                         randomInit=randomInit,
                         maxNumOptIter=maxNumOptIter,
                         numOptRelTol=numOptRelTol, init=init, method=method, hessian=hessian,
                         learning.rate=learning.rate)
      o <- order(ans$Mixture$Weight, decreasing=TRUE)
      ans <- within(ans, {
          Group <- Group[,o, drop=FALSE]
          Mixture$Weight <- Mixture$Weight[o]
          Fit <- rapply(Fit, function(mat){mat[, o, drop=F]}, how='replace')
          EM.diagnostics <- EM.diagnostics
          Data <- Data
          nll.data <- nll.data
          nLB.data <- nLB.data
          Ez <- Ez
          EM_lambda_optim_message <- EM_lambda_optim_message
          EM_alpha_list <- EM_alpha_list
          Ez_list <- Ez_list
      })
      with(ans, .DMN(goodnessOfFit=GoodnessOfFit,
                     group=Group,
                     mixture=Mixture,
                     fit=Fit,
                     EM.diagnostics=EM.diagnostics,
                     Data=Data,
                     nll.data=nll.data,
                     nLB.data=nLB.data,
                     Ez=Ez,
                     EM_lambda_optim_message=EM_lambda_optim_message,
                     EM_alpha_list=EM_alpha_list,
                     Ez_list=Ez_list))
    } #repet.func ends

    if (length(K) == 1) {
      if (parallel && (repetition > 1)) {
        results <- mclapply(seq_len(repetition),
                            repet.func,
                            k=K,
                            mc.cores=getOption("mc.cores", min(repetition, 4L)))
      } else if (parallel){
        results <- lapply(seq_len(repetition), repet.func, k=K, parallel.shift=T)
      } else {
        results <- lapply(seq_len(repetition), repet.func, k=K)
      }

      best <- which.min(sapply(results, BIC))
      return(results[[best]])
    } #length(K) == 1
    else {
      k.and.r <- expand.grid(K=K, R=seq_len(repetition))
      if (parallel) {
        results <- mclapply(seq_len(nrow(k.and.r)),
                          function(row.i){
                            k <- k.and.r[row.i, 'K']
                            r <- k.and.r[row.i, 'R']
                            repet.func(r, k)
                          },
                          mc.cores=min(length(K)*repetition, detectCores()))
      } else {
        results <- lapply(seq_len(nrow(k.and.r)),
                          function(row.i){
                            k <- k.and.r[row.i, 'K']
                            r <- k.and.r[row.i, 'R']
                            repet.func(r, k)
                            })
      }
      
      results <- split(results, k.and.r$K)
      results <- lapply(results, function(rep)rep[[which.min(sapply(rep, BIC))]])
      return (results)
    }
}

mixture <-
    function(object, ..., assign=FALSE)
{
    if (assign) {
        apply(mixture(object), 1, which.max)
    } else {
        object@group
    }
}

## Dirichlet

goodnessOfFit <- function(object, ...) object@goodnessOfFit

laplace <- function(object, ...) goodnessOfFit(object)[["Laplace"]]

.AIC.DMN <- function(object, ...) goodnessOfFit(object)[["AIC"]]

setMethod(AIC, "DMN", .AIC.DMN)

.BIC.DMN <- function(object, ...) goodnessOfFit(object)[["BIC"]]

setMethod(BIC, "DMN", .BIC.DMN)

mixturewt <-
    function(object, ...)
{
    if (is.list(fitted(object)))
      list(pi=object@mixture$Weight, theta=lapply(fitted(object), colSums))
    else  #for backward compatibility
      data.frame(pi=object@mixture$Weight, theta=colSums(fitted(object)))
}

.fitted.DMN <- function(object, ..., scale=FALSE)
{
    fit <- object@fit$Estimate
    if (scale) {
      if (is.list(fit))
        fit <- lapply(fit, function(f) scale(f, FALSE, colSums(f)))
      else #backward compatibility
        fit <- scale(fit, FALSE, mixturewt(object)$theta)
    }
    fit
}

setMethod(fitted, "DMN", .fitted.DMN)

## predict

.neg_log_evidence_i <-
    function(x, alpha)
{
    .B <- function(x)
        sum(lgamma(x)) - lgamma(sum(x))
    -(.B(x + alpha) - .B(alpha))
}

.predict.DMN <-
    function(object, newdata, ..., logevidence=FALSE)
{
    if (is.vector(newdata))
        newdata <- matrix(newdata, nrow=1)
    lambda <- fitted(object)

    K <- ncol(lambda)
    alpha <- sapply(seq_len(K), function(k, lamda, x) {
        apply(x, 1, .neg_log_evidence_i, lambda[,k])
    }, lambda, newdata)
    if (is.vector(alpha))
        alpha <- matrix(alpha, nrow=1,
                        dimnames=list(rownames(newdata), NULL))
    if (!logevidence) {
        wt <- mixturewt(object)$pi
        offset <- apply(alpha, 1, min)
        z <- sweep(exp(-(alpha - offset)), 2, wt, "*")
        z / rowSums(z)
    } else {
        alpha
    }
}

setMethod(predict, "DMN", .predict.DMN)

## print / plot

setMethod(show, "DMN",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("K:", ncol(mixture(object)), "\n")
    if (is.list(fitted(object))) {
      cat("M:", length(fitted(object)), "\n")
      cat("rows x cols: ", paste(nrow(mixture(object)), sapply(fitted(object), nrow), sep='x'), "\n")
    } else {
      cat("rows x cols: ", nrow(mixture(object)), "x", nrow(fitted(object)), "\n")
    }
    cat("Laplace:", laplace(object), "BIC:", BIC(object),
       "AIC:", AIC(object), "\n")
})
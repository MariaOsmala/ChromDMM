##  utility

.gradient <-                   # RColorBrewer::brewer.pal(9, "YlOrRd")
    c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
      "#FC4E2A", "#E31A1C", "#BD0026", "#800026")

.divergent <-                  # RColorBrewer::brewer.pal(9, "RdYlBu")
    c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
      "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")

.qualitative <-                # RColorBrewer::brewer.pal(10, "Paired")
    c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

csubset <- function(val, x, pheno, cidx=TRUE)
{
    ridx <- pheno %in% val
    if (!cidx)
        cidx <- colSums(x[ridx,]) != 0
    x[ridx, cidx]
}

#called from c++ shift_and_flip_signal(binned=binned.data,
#                      indices=1
#                      dist=dist.candidates[1],
#                      flip=F;
#
#this function shifts the samples (given as indices) by a shift=dist
#this function is also used to flip the samples
shift.and.flip.signal <- function(binned, indices, dist, flip) {
  stopifnot(length(indices) == length(dist))
  stopifnot(length(indices) == length(flip))

  if (length(indices) == 0)
    return(binned)

  #shift IRangesList windows
  windows <- attr(binned, 'windows') #list of IRanges objects 121-1880 (start-end) width 1760
  indseq <- rep(0, length(windows))
  indseq[indices] <- -dist #shift window backwards to shift reads forward
  windows <- shift(windows, indseq) #only the indices(th) window shifted
  attr(binned, 'windows') <- windows

  #update shifts vector
  shifts <- attr(binned, 'shifts')
  shifts[indices] <- shifts[indices] + dist
  attr(binned, 'shifts') <- shifts

  #update limits, the possible limits change to prevent exceeding the boundaries
  limits <- attr(binned, 'shift.limits')
  limits[indices,] <- limits[indices,] - dist
  attr(binned, 'shift.limits') <- limits

  #update flips
  flips <- attr(binned, 'flips')
  flips[indices] <- flip
  attr(binned, 'flips') <- flips

  #finally extract new shifted data
  binned <- extract_binned_signal(binned, indices)

  binned
}

#' Most sophisticated integer minimizer
#'
#' @param fn optimization_func,
#' @param limits 
#' @param index sample indexes?
#' @param step 
#' @param parallel 
#' @param verbose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
most.sophisticated.integer.minimizer <- function(fn,
                                                 limits,
                                                 index,
                                                 step=10,
                                                 parallel=F,
                                                 verbose=F,
                                                 ...) {
  initparams <- rep(0, nrow(limits))

  fminimizer <- function(pindex, pb){
    if(!is.null(pb)) setTxtProgressBar(pb, pindex)

    minlimit <- limits[pindex, 1]
    maxlimit <- limits[pindex, 2]
    p <- initparams[pindex]
    i <- index[pindex]

    if (((p-step) < minlimit) || ((p+step) > maxlimit))
      return(p)

    negresult <- fn(p-step, i, ...)
    posresult <- fn(p+step, i, ...)
    initialresult <- fn(p, i, ...)
    
    #optimization_func(IntegerVector shift_dist,
    #                  IntegerVector index,
    #                  LogicalVector flip,
    #                  List data,
    #                  List alpha,
    #                  NumericMatrix Z)
    #negresult <- fn(p-step, i, flip=F,data, alpha, Z)
    #posresult <- fn(p+step, i, flip=F,data, alpha, Z)
    #initialresult <- fn(p, i, flip=F,data,alpha,Z)

    if (initialresult < negresult && initialresult < posresult) #shifting does not help
      return(p)

    if (negresult < posresult) {
      progfn <- `-`
      result <- negresult
    } else {
      progfn <- `+`
      result <- posresult
    }
    prevresult <- result
    p <- progfn(p, step)

    repeat {

      nextp <- progfn(p, step)

      if (nextp < minlimit || nextp > maxlimit) #we are over limits
        break

#       print(sprintf('Trying %d...', nextp))
      result <- fn(nextp, i, ...)
      #result <- fn(nextp, i, flip, data, alpha,Z)
      if(prevresult <= result) #we have reached the optimum
        break

      prevresult <- result
      p <- nextp
    } #repeat
    p
  } #fminimizer
  
  if(parallel) {
    simplify2array(mclapply(seq_along(initparams), fminimizer, pb=NULL, mc.cores=detectCores()))
  } else {
    if (verbose)
      pb <- txtProgressBar(1, length(initparams),
                           initial=1, title='Shift optimization',
                           width=40,
                           style=3)
    else
      pb <- NULL

    sapply(seq_along(initparams), fminimizer, pb=pb)
  }
}


#' brute.force.integer.minimizer
#'
#' try all possible shifts and flips and minimize given optimization function
#'
#' @param fn optimization_func, what do we optimize???
#' @param limits attr(binned.data, 'shift.limits'), N x 2 (min and max)
#' @param step =bin.width (40)
#' @param parallel.shift (default F)
#' @param verbose
#' @param shift.reads T or F
#' @param flip T or F
#' @data  =binned.data,
#' @alpha =lapply(lambda, exp)
#' @Z     =Ez
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'

# fn=optimization_func
# limits=attr(binned.data, 'shift.limits')
# step=bin.width
# parallel.shift=parallel.shift
# verbose=verbose
# shift.reads=shift.reads
# flip=flip
# data=binned.data
# alpha=lapply(lambda, exp)
# Z=Ez


brute.force.integer.minimizer <- function(fn,
                                          limits,
                                          step=10,
                                          parallel.shift=F,
                                          verbose=F,
                                          shift.reads=T,
                                          flip=F,
                                          ...) {
  initparams <- rep(0, nrow(limits)) #N

  stopifnot(shift.reads || flip)
  #res <- lapply(seq_along(initparams), fminimizer, progress.bar=progress.bar)
  #pindex is 1:N
  #fminimizer(1, progress.bar)
  fminimizer <- function(pindex, progress.bar){
    if(!is.null(progress.bar)) setTxtProgressBar(progress.bar, pindex)

    if(shift.reads) {
      minlimit <- limits[pindex, 1]
      maxlimit <- limits[pindex, 2]

      #dist.candidates <- seq(minlimit, maxlimit, step)
      #-120  -80  -40    0   40   80  120
      dist.candidates <- c(rev(seq(0, minlimit, -step)), seq(0, maxlimit, step)[-1])
    }
    else {
      dist.candidates <- 0
    }

    best <- list(shift=0, flip=F)

    #first shifts
    #optimization_func( dist.candidates[1], pindex, flip=F)
    result <- sapply(dist.candidates, function(p)fn(p, pindex, flip=F, ...))
    
    # result <- sapply(dist.candidates, function(p)fn(p, pindex, flip=F, data, alpha, Z))
    if(flip) {
      result.with.flip <- sapply(dist.candidates, function(p)fn(p, pindex, flip=T, ...))
      #result.with.flip <- sapply(dist.candidates, function(p)fn(p, pindex, flip=T, data, alpha, Z))
      if (min(result.with.flip) < min(result)) { #flipping improves the value
        best$flip <- T
        best$shift <- dist.candidates[which.min(result.with.flip)]
      } else { #flipping does not improve the value
        best$shift <- dist.candidates[which.min(result)]
        best$flip <- F
      }
    } else {
      best$shift <- dist.candidates[which.min(result)]
      best$flip <- F
    }

    #print(sprintf('Best distance for %d is %d...', pindex, best))

    return(best)
  }

  if(parallel.shift) {
    res <- mclapply(seq_along(initparams),
                    fminimizer,
                    progress.bar=NULL,
                    mc.cores=detectCores())
  } else {
    if (verbose) {
      progress.bar <- txtProgressBar(1, length(initparams),
                           initial=1, title='Shift/flip optimization',
                           width=40,
                           style=3)
    } else {
      progress.bar <- NULL
    }

    res <- lapply(seq_along(initparams), fminimizer, progress.bar=progress.bar)
  }

  return(list(shift=sapply(res, function(x)x$shift),
              flip=sapply(res, function(x)x$flip)))
}


# shift.optimization.func <- function(d, index, data, alpha, Z) {
#   library(gsl)
#   shifteddata <- shift.signal(data, index, d)
#
#   M <- length(shiftedata)
#   K <- nrow(Z)
#   ret <- 0
#
#   for (k in seq_len(K)) {
#     for (m in seq_len(M)) {
#       datasub <- shifteddata[[m]][index,]
#       xsumalpha <- t(t(datasub) + alpha[[m]][k,])
#       betaXalpha <- rowSums(lngamma(xsumalpha)) - lngamma(rowSums(xsumalpha))
#       logsumX <- lnfact(rowSums(datasub))
#       factlogsum <- rowSums(lnfact(datasub))
#       ret <- ret + sum(Z[k, index]*(betaXalpha + logsumX - factlogsum))
#     }
#   }
#   -ret
# }

# same function without using gsl package
# shift.optimization.func2 <- function(d, index, data, alpha, Z) {
#
#   shifteddata <- shift.signal(data, index, d)
#
#   M <- length(data)
#   K <- nrow(Z)
#   ret <- 0
#
#   for (k in seq_len(K)) {
#     for (m in seq_len(M)) {
#       datasub <- data[[m]][index,]
#       xsumalpha <- t(t(datasub) + alpha[[m]][k,])
#       betaXalpha <- rowSums(lgamma(xsumalpha)) - lgamma(rowSums(xsumalpha))
#       logsumX <- lfactorial(rowSums(datasub))
#       factlogsum <- rowSums(lfactorial(datasub))
#       ret <- ret + sum(Z[k, index]*(betaXalpha + logsumX - factlogsum))
#     }
#   }
#   -ret
# }
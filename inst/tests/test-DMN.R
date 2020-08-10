library(testthat)
library_if_available(DirichletMultinomial)

context("Testing shifting optimization function()")

N <- 2
K <- 2
M <- 2

X <- list(FirstDatatype=rbind(c(0,2,0),
                              c(1,1,1)),
          SecondDatatype=rbind(c(0,0,1),
                               c(1,7,1)))

alpha <- list(rbind(c(7.56, 8.49, 1.12),
                    c(0.02, 7.1, 0.004)),
              rbind(c(1,2,3),
                    c(0.5, 1.5, 2.5)))

Z <- rbind(c(0.01, 0.6),
           c(0.99, 0.4))

logbeta <- function(vec) {sum(lgamma(vec))-lgamma(sum(vec))}

## test no permutation
test_that("Shift optimization func working correctly", { #ikm
  exp <- -sum(0.01*(sum(logbeta(X[[1]][1,]+alpha[[1]][1,])) + log(factorial(sum(X[[1]][1,]))) - sum(log(factorial(X[[1]][1,])))), #111
            0.99*(sum(logbeta(X[[1]][1,]+alpha[[1]][2,])) + log(factorial(sum(X[[1]][1,]))) - sum(log(factorial(X[[1]][1,])))), #121
            0.01*(sum(logbeta(X[[2]][1,]+alpha[[2]][1,])) + log(factorial(sum(X[[2]][1,]))) - sum(log(factorial(X[[2]][1,])))), #112
            0.99*(sum(logbeta(X[[2]][1,]+alpha[[2]][2,])) + log(factorial(sum(X[[2]][1,]))) - sum(log(factorial(X[[2]][1,])))), #122
            0.6*(sum(logbeta(X[[1]][2,]+alpha[[1]][1,])) + log(factorial(sum(X[[1]][2,]))) - sum(log(factorial(X[[1]][2,])))), #211
            0.4*(sum(logbeta(X[[1]][2,]+alpha[[1]][2,])) + log(factorial(sum(X[[1]][2,]))) - sum(log(factorial(X[[1]][2,])))), #221
            0.6*(sum(logbeta(X[[2]][2,]+alpha[[2]][1,])) + log(factorial(sum(X[[2]][2,]))) - sum(log(factorial(X[[2]][2,])))),#212
            0.4*(sum(logbeta(X[[2]][2,]+alpha[[2]][2,])) + log(factorial(sum(X[[2]][2,])))  - sum(log(factorial(X[[2]][2,]))))) #222

  expect_that(optimization_func(d=c(0,0),
                                data = X,
                                alpha = alpha,
                                Z = Z), equals(exp))
})
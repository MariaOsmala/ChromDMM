library("ChromDMM")
library("optparse")

# options
option_list = list(
  make_option(c("-d", "--data"), action="store", default=NULL, type='character',
              help="The file containing the data to cluster."),
  make_option(c("-c", "--cluster"), action="store", default="1", type='character',
              help="The number of clusters. Multiple cluster numbers can be given, separated by comma, e.g. 1,2,3"),
  make_option(c("--bin.size", "-B"), action="store", default=1, type='numeric',
              help="Bin size if the data resolution needs to be decreased. No change in bin size in default."),
  make_option(c("--verbose"), action="store_true", default=TRUE,
              help="Print the progress of the EM algorithm. Set to FALSE, if --parallel==TRUE"),
  make_option(c("--shift", action="store", default=1, type="numeric"),
              help="Number of shift states, by default 1 (corresponds to no shifting)."),
  make_option(c("--flip"), action="store_true", default=FALSE,
              help="Enables flipping"),
  make_option(c("--seed.boolean"), action="store_true", default=FALSE, 
              help="Whether to use seed for the random number generator"),
  make_option(c("--seed"), action="store", default=NULL, type='integer',
              help="An integer to set the random number generator seed. Used only if seed-boolean==TRUE"),
  make_option(c("--eta"), action="store", default=1.1, type='numeric',
              help="Dirichlet parameter Gamma hyperprior shape"),
  make_option(c("--nu"), action="store", default=0.1, type='numeric',
              help="Dirichlet parameter Gamma hyperprior rate."),
  make_option(c("--etah"), action="store", default=10, type='numeric',
              help="Regularisaion parameter Gamma hyperprior shape."),
  make_option(c("--nuh"), action="store", default=10, type='numeric',
              help="Regularisation parameter Gamma hyperprior rate"),
  make_option(c("--EM.max.iter"), action="store", default=250, type='numeric',
              help="The maximal number of EM iterations."),
  make_option(c("--EM.num.tol"), action="store", default=1e-6, type='numeric',
              help="The numerical tolerance to check the convergence of EM algorithm"),
  make_option(c("--BFGS.max.iter"), action="store", default=1000, type='numeric',
              help="The maximal number of BFGS iterations."),
  make_option(c("--BFGS.num.tol"), action="store", default=1e-6, type='numeric',
              help="The numerical tolerance to check the convergence of BFGS"),
  make_option(c("--initialisation"), action="store", default="kmeans++", type='character',
              help="The initialisation strategy for the EM algorith. Options: kmeans++ (default) and random."),
  make_option(c("--soft.kmeans.maxit"), action="store", default=1000, type='numeric',
              help="The maximal number of soft-kmeans iterations in the EM initialisation."),
  make_option(c("--soft.kmeans.stiffness"), action="store", default=50, type='numeric',
              help="Parameter for the soft kmeans in the EM initialisation"),
  make_option(c("--soft.kmeans.randomInit"), action="store_true", default=TRUE, 
              help="Is the soft kmeans initialised randomly (default TRUE)."),
  make_option(c("--repetition"), action="store", default=1, type='numeric',
              help="The number of EM algorithm repetitions. Each repetition is initialised randomly (seed-boolean needs to be FALSE), 
              and the best fit of the repetitions (according to BIC) is chosen as the final model fit."),
  make_option(c("--parallel"), action="store_true", default=TRUE, 
              help="Whether the EM optimisations for varying number of clusters and for the repetitions are performed in parallel."),
  make_option(c("--output"), action="store", default="", type='character',
              help="The file where the results will be saved (in RDS format).")
)





opt = parse_args(OptionParser(option_list=option_list))

opt$shift=as.numeric(opt$shift)

if(opt$output == "")
{ stop("Error! Output file not specified (--output)!") }

if(opt$data == "")
{ stop("Error! Data file not specified (--data)!") }

if(opt$shift%%2 == 0)
{ stop("Error! Number of shift states needs to be an odd number (--shift)!") }

if(opt$parallel == TRUE && opt$verbose==TRUE)
{ print("Note! if parallel=TRUE, verbose is set to FALSE") 
  opt$verbose==FALSE
  }

# print(paste0("data: ", str(opt$data)))
# opt$cluster=as.numeric(unlist(strsplit(opt$cluster, ",")))
# print(paste0("cluster: ", str(opt$cluster)))
# print(paste0("bin.size: ",str(opt$bin.size)))
# print(paste0("verbose: ",str(opt$verbose)))
# print(paste0("shift: ", str(opt$shift)))
# print(paste0("flip: ", str(opt$flip)))
# print(paste0("seed-boolean: ", str(opt$seed.boolean)))
# print(paste0("seed: ", str(opt$seed)))
# print(paste0("eta: ", str(opt$eta)))
# print(paste0("nu: ", str(opt$nu)))
# print(paste0("etah: ", str(opt$etah)))
# print(paste0("nuh: ", str(opt$nuh)))
# print(paste0("EM.max.iter: ",str(opt$EM.max.iter)))
# print(paste0("EM.num.tol: ", str(opt$EM.num.tol)))
# print(paste0("BFGS.max.iter: ",str(opt$BFGS.max.iter)))
# print(paste0("BFGS.num.tol: ",str(opt$BFGS.num.tol)))
# print(paste0("initialisation: ",str(opt$initialisation)))
# print(paste0("soft.kmeans.maxit: ",str(opt$soft.kmeans.maxit)))
# print(paste0("soft.kmeans.stiffness: ",str(opt$soft.kmeans.stiffness)))
# print(paste0("soft.kmeans.randomInit: ",str(opt$soft.kmeans.randomInit)))
# print(paste0("repetition: ",str(opt$repetition)))
# print(paste0("parallel: ",str(opt$parallel)))
# print(paste0("output: ", str(opt$output)))

print(paste0("data: ", opt$data))
opt$cluster=as.integer(unlist(strsplit(opt$cluster, ",")))
print(paste0("cluster: ", opt$cluster))
print(str(opt$cluster))
print(str(1:3))
print(paste0("bin.size: ",opt$bin.size))
print(paste0("window: ", opt$window))
print(paste0("verbose: ",opt$verbose))
print(paste0("shift: ", opt$shift))
print(paste0("flip: ", opt$flip))
print(paste0("seed-boolean: ", opt$seed.boolean))
print(paste0("seed: ", opt$seed))
print(paste0("eta: ", opt$eta))
print(paste0("nu: ", opt$nu))
print(paste0("etah: ", opt$etah))
print(paste0("nuh: ", opt$nuh))
print(paste0("EM.max.iter: ",opt$EM.max.iter))
print(paste0("EM.num.tol: ", opt$EM.num.tol))
print(paste0("BFGS.max.iter: ",opt$BFGS.max.iter))
print(paste0("BFGS.num.tol: ",opt$BFGS.num.tol))
print(paste0("initialisation: ",opt$initialisation))
print(paste0("soft.kmeans.maxit: ",opt$soft.kmeans.maxit))
print(paste0("soft.kmeans.stiffness: ",opt$soft.kmeans.stiffness))
print(paste0("soft.kmeans.randomInit: ",opt$soft.kmeans.randomInit))
print(paste0("repetition: ",opt$repetition))
print(paste0("parallel: ",opt$parallel))
print(paste0("output: ", opt$output))

# 




# opt=list()
# opt$data="experiment_data/data.RDS"
# opt$cluster=2
# opt$bin.size=1
# opt$window=50
# opt$verbose=TRUE
# opt$shift=21
# opt$flip=TRUE
# opt$seed.boolean=FALSE
# opt$seed=FALSE
# opt$eta=1.1
# opt$nu=0.1
# opt$etah=10
# opt$nuh=10
# opt$EM.max.iter=250
# opt$EM.num.tol=1e-06
# opt$BFGS.max.iter=1000
# opt$BFGS.num.tol=1e-06
# opt$initialisation="kmeans++"
# opt$soft.kmeans.maxit=1000
# opt$soft.kmeans.stiffness=50
# opt$soft.kmeans.randomInit=TRUE
# opt$repetition=1
# opt$parallel=TRUE
# opt$output="data_experiments/results.RDS"




seed=FALSE
if(opt$seed.boolean==TRUE){
  seed=opt$seed
}


#Artificial data
zeta=NULL
flip=FALSE

if(opt$flip==TRUE){
  zeta=t(as.matrix(c(0.5,0.5)))
  flip=TRUE
}

xi=NULL
shift.reads=FALSE

if(opt$shift>1){ #check also that opt$shift is an odd number
  shift.reads=TRUE
  pyramid <- function(S){
    xi=rep(0,S)
    xi[floor(S/2)+1]=floor(S/2)+1
    for(i in floor(S/2):1){
      xi[i]=i
      xi[S-i+1]=i
    }
    xi=xi/sum(xi)
  }
  
  xi=pyramid(opt$shift)
  
  xi=t(as.matrix(xi))
  
}


data<- readRDS( opt$data) 

#Remove samples which are zero for some modification? almost everything get's removed?
zero.indices <- which(apply(sapply(data$data, rowSums), 
                            1, prod) == 0)

if (length(zero.indices) > 0) {
  warning(sprintf('Excluding are all-zero row(s) in the data: %s\n',
                  paste(names(zero.indices), collapse = ',')))
}

fit <- dmn(count=data$data,
    K=opt$cluster,
    bin.width=opt$bin.size,
    S=opt$shift,
    verbose=opt$verbose,
    seed=seed,
    shift.reads=shift.reads,
    flip=opt$flip,
    eta=opt$eta,
    nu=opt$nu,
    etah=opt$etah,
    nuh=opt$nuh,
    maxIt=opt$EM.max.iter,
    xi=xi,
    zeta=zeta,
    EM.threshold=opt$EM.num.tol,
    soft.kmeans.maxit=opt$soft.kmeans.maxit,
    soft.kmeans.stiffness=opt$soft.kmeans.stiffness,
    randomInit=opt$soft.kmeans.randomInit,
    repetition=opt$repetition,
    maxNumOptIter=opt$BFGS.max.iter,
    numOptRelTol=opt$BFGS.num.tol,
    parallel=opt$parallel, init=opt$initialisation, hessian=FALSE)

# maxNumOptIter=opt$BFGS.max.iter
# numOptRelTol=opt$BFGS.num.tol
# maxIt=opt$EM.max.iter
# soft.kmeans.maxit=opt$soft.kmeans.maxit
# soft.kmeans.stiffness=opt$soft.kmeans.stiffness
# randomInit=opt$soft.kmeans.randomInit
# repetition=opt$repetition
# parallel=opt$parallel
# init=opt$initialisation
# EM.threshold=opt$EM.num.tol




saveRDS(fit, opt$output)





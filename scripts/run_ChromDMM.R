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
  make_option(c("--window", "-W"), action="store", default=1, type='numeric',
              help="Window size if adjusting needed. No change if the window size matches with the data window size."),
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

# print(paste0("data: ", str(opt$data)))
# opt$cluster=as.numeric(unlist(strsplit(opt$cluster, ",")))
# print(paste0("cluster: ", str(opt$cluster)))
# print(paste0("bin.size: ",str(opt$bin.size)))
# print(paste0("window: ", str(opt$window)))
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
opt$cluster=as.numeric(unlist(strsplit(opt$cluster, ",")))
print(paste0("cluster: ", opt$cluster))
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



print("set seed")
seed=FALSE
if(opt$seed.boolean==TRUE){
  seed=opt$seed
}


#Artificial data
zeta=NULL
flip=FALSE
print("flip")
if(opt$flip==TRUE){
  zeta=t(as.matrix(c(0.5,0.5)))
  flip=TRUE
}

xi=NULL
shift.reads=FALSE
print("shift")
print(str(opt$shift))
if(opt$shift>1){ #check also that opt$shift is an odd number
  shift.reads=TRUE
  print("shift2")
  pyramid <- function(S){
    xi=rep(0,S)
    xi[floor(S/2)+1]=floor(S/2)+1
    for(i in floor(S/2):1){
      xi[i]=i
      xi[S-i+1]=i
    }
    xi=xi/sum(xi)
  }
  print("shift3")
  xi=pyramid(opt$shift)
  print("shift4")
  xi=t(as.matrix(xi))
  
}

print("read data")
data<- readRDS( opt$data) 
print("fit model")
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

stop()

#TRUE data
################## Enhancers 10 modifications ##########################


#41) ARGS="1000_enhancers_bin_1_window_2000_only5prime.RData 
# 40 
# dmn-result-1000-enhancers-withShift-withFlip_40_10mods_withoutRep_1.Rds
# enhancers
# 10 
# H2AZ 
# H3K27ac 
# H3K4me1 
# H3K4me2 
# H3K4me3 
# H3K79me2 
# H3K9ac 
# RNAPOL2 DNase-seq MNase-seq 8 1.1 0.1 10 10 TRUE TRUE 1 TRUE FALSE 1000";; #

args = commandArgs(trailingOnly=TRUE)
signals_name<-as.character(args[1])
bin_size <-as.numeric(args[2])
save_name<-as.character(args[3])
region_type=as.character(args[4])
mod_number=as.numeric(args[5])
mod_subset=c()
j=6
for(i in 1:mod_number){
  mod_subset[i]=as.character(args[5+i])
  j=j+1
}

Kmax=as.numeric(args[j]) #The max number fo Dirichlet components to fit
eta=as.numeric(args[j + 1])
nu=as.numeric(args[j + 2])
etah=as.numeric(args[j + 3])
nuh=as.numeric(args[j + 4])
shift.reads=as.logical(args[j + 5])
flip=as.logical(args[j + 6])
repetition=as.numeric(args[j + 7])
parallel=as.logical(args[j + 8])

directed=as.logical(args[j + 9])
maxIt=as.numeric(args[j + 10])

N=as.numeric(args[j + 11])

print(signals_name)
print(bin_size)
print(save_name)
print(region_type)
print(mod_number)
print(mod_subset)
print(Kmax)
print(eta)
print(nu)
print(etah)
print(nuh)
print(shift.reads)
print(flip)
print(repetition)
print(parallel)
print(directed)
print(maxIt)

save_name=gsub(".Rds", paste0("_EMiter_",maxIt,".Rds"), save_name)
#if signals created by Göcken, signals.CTCF.enhancer.signal or promoter.signal

#load(paste0("/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/enhancerPrediction2020/Data/K562/",signals_name))
load(paste0("/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/preprint_devel/Data/K562/data_R/",signals_name))
#load(paste0("/Volumes/scratch/cs/csb/projects/enhancer_prediction/experiments/RProjects/preprint_devel/Data/K562/data_R/",signals_name))
if(region_type=="enhancers"){
  signals=normalized_profiles  
}
if(region_type=="promoters"){
  if(directed==TRUE){
    signals=normalized_profiles_directed  
  }else{
    signals=normalized_profiles_undirected  
  }
  
}

for(s in names(signals)){
  signals[[s]]=signals[[s]][,1:N]
  
}


non_histone_ind=c(grep( 'Dnase',names(signals)),grep( 'Nsome',names(signals)))
clearer_names<-names(signals)
clearer_names=unlist(sapply(clearer_names,strsplit, "RawData"))
clearer_names[non_histone_ind]=c("DNase-seq", "MNase-seq")
clearer_names[non_histone_ind]=c("DNase-seq","MNase-seq")
clearer_names=sub("Ctcf","CTCF",clearer_names)
clearer_names=sub("az","AZ",clearer_names)
clearer_names=sub("Pol2","RNAPOL2",clearer_names)
clearer_names=sub("k","K",clearer_names)
names(clearer_names)=NULL
names(signals)=clearer_names



save_path="RData/"
figure_path="figures/"

print(paste0(save_path, save_name))
print(getwd())
library("RcppGSL")
library("DMMChrom")

# Bin signals ============================================================

signals.binned <- Map( function(mod, modname) {
  count <- bin_signal(t(as.matrix(mod)), 1)
  #str(count)
  rownames(count) <- paste0(modname, '.', seq(nrow(count)))
  colnames(count) <- paste0('bin', seq(ncol(count)))
  
  count[which(count<0)]=0
  count=round(count)
  
  
  count
  
}, signals, names(signals))


signals.binned2 <- Map( function(mod, modname) {
  count <- bin_signal(t(as.matrix(mod)), bin_size)
  #str(count)
  rownames(count) <- paste0(modname, '.', seq(nrow(count)))
  colnames(count) <- paste0('bin', seq(ncol(count)))
  
  count[which(count<0)]=0
  count=round(count)
  
  
  count
  
}, signals, names(signals))



# enhancers.binned <- Map(function(mod, modname) {
#   count <- bin_signal(t(as.matrix(mod$enhancer.signal)), 40)
#   rownames(count) <- paste0(modname, '.', seq(nrow(count)))
#   colnames(count) <- paste0('bin', seq(ncol(count)))
#   count
# 
# }, signals, names(signals))

signals.subset <- signals.binned[mod_subset]
signals.subset2<-signals.binned2[mod_subset]


#Study the distribution

#library("ggplot2")
#p<-gList()
#gamma_param=data.frame(shape=rep(0,length(mod_subset)), rate=rep(0,length(mod_subset)),row.names=mod_subset)
#for(mod in mod_subset){
#  df<- data.frame(colMean=rowMeans(signals.subset2[[mod]]))
#  gamma_param[mod,]<-fitdistr(df$colMean, densfun="gamma")$estimate
#  p[[mod]] <- ggplotGrob(ggplot(df, aes(x=colMean))+geom_histogram(binwidth=1)
#                        +ggtitle(mod)+xlab("Mean")+ylab("Frequency"))
#
#  print(p[[mod]])
#}
#library("gridExtra")
#png("test.png",
#    width=500,
#    height=700)
#marrangeGrob(p, nrow=5, ncol=2)
#dev.off()

#Fit gamma distribution into these values



# Parallel ============================================================
#signals.subset is a list of matrixes of size 1000 x 50

# count=signals.subset
# K=1:Kmax 
# bin.width = bin_size
# shift.ratio = 1/8
# verbose = FALSE
# seed = F
#shift.reads = T
#flip = F 
#eta = 0.1
#nu = 0.1, 
#etah = 0
#nuh = 0
# maxIt = 250
# EM.threshold = 1e-06
# soft.kmeans.maxit = 1000
# soft.kmeans.stiffness = 50
# randomInit = T
# maxNumOptIter = 1000
# numOptRelTol = 1e-12
#parallel = T

#dmn bins the data signals.subset should be in 1 bp format

#Remove samples which are zero for some modification? almost everything get's removed?
zero.indices <- which(apply(sapply(signals.subset, rowSums), 
                            1, prod) == 0)

if (length(zero.indices) > 0) {
  warning(sprintf('Excluding are all-zero row(s) in the data: %s\n',
                  paste(names(zero.indices), collapse = ',')))
}

#save.image("/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/DMM-private-master-devel-works/startingPoint.RData")

#save in Spar-K compatible format(needs counts, what about normalized?)
#for(mod in mod_subset){
#write.table(signals.subset2[[mod]],quote=FALSE, sep=" ", row.names=FALSE, 
#            col.names=FALSE,
#            file=paste0("/m/cs/scratch/csb/projects/enhancer_clustering/EnhancerClustering/method_comparison/Spar-K/",mod,".txt"))
#}


zeta=NULL
S=1
xi=NULL
if(flip==TRUE){
  zeta=t(as.matrix(c(0.5,0.5)))
}
if(shift.reads==TRUE){
  S=21
  pyramid <- function (S){
    xi=rep(0,S)
    xi[floor(S/2)+1]=floor(S/2)+1
    for(i in floor(S/2):1){
      xi[i]=i
      xi[S-i+1]=i
    }
    xi=xi/sum(xi)
  }
  xi=pyramid(S)
  xi=t(as.matrix(xi))
}


try.catch=try(rnorm(1))
class(try.catch)="try-error"

#concatenate and write the data into Spar-K and ChIP-Partitioning compatible format
create_spark_data=FALSE
if(create_spark_data==TRUE){
  spark.data <- do.call(cbind, signals.subset)
  attr(spark.data, "dimnames")[[2]]= paste0('bin', seq_len(ncol(spark.data)))
  attr(spark.data, "dimnames")[[1]]= paste0('loci', seq_len(nrow(spark.data)))
  write.table(spark.data, file=paste0("real_enhancer_data_spark_chippart.txt"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
}


while(class(try.catch)=="try-error"){
  try.catch <- try( subset.fits <- ChromDMM::dmn(count=signals.subset, K=1:Kmax, 
                                                 bin.width=bin_size, 
                                                 verbose = F, eta=eta, nu=nu, etah=etah, nuh=nuh, S=S,
                                                 shift.reads = shift.reads, flip=flip, zeta=zeta, xi=xi,
                                                 repetition=repetition, randomInit=T, parallel=parallel,EM.threshold=1e-6, 
                                                 maxNumOptIter = 1000, numOptRelTol=1e-6, maxIt=maxIt, init="kmeans++" )
  )
}
print(getwd())
print(paste0(save_path, save_name))
saveRDS(subset.fits, paste0(save_path, save_name))
saveRDS(zero.indices, paste0(save_path, gsub(".Rds", "_zeroindices.Rds", save_name)))



#  opt=list()
#  data_path="/Users/mpirttin/projects/EnhancerClustering/experiments/MariasExperimentsPackageProperlyInstalled/enhancerdnase-experiments/artificial_data_better/chippartitioning/"
# # 
# # # 
#  n_first=10
#  n_second=10
# # # 
# # #for reference
# # opt$data=paste0(data_path,"shift-concat-2-2-",n_first,"-",n_second,"-eta-1.1-nu-0.1-etah-1.1-nuh-0.1-not-shifted-bin-40.txt")
# # 
#  opt$data=paste0(data_path,"shift-concat-2-2-",n_first,"-",n_second,"-eta-1.1-nu-0.1-etah-1.1-nuh-0.1-shifted-bin-40.txt")
# #               
# opt$cluster=2
# opt$iter=100
# opt$shift=21
# opt$seed=1
# opt$output=paste0("results_artificial_data_better/shift-concat-2-2-",n_first,"-",n_second,"-eta-1.1-nu-0.1-etah-1.1-nuh-0.1/result.Rds")


# parses options
opt = parse_args(OptionParser(option_list=option_list))
opt$shift = as.numeric(opt$shift)
if(opt$output == "")
{ stop("Error! Output file not specified (--output)!") }


# read data
data = as.matrix(read.table(opt$data))

print(opt$seeding)
# sets the random number generator
set.seed(opt$seed)

if(opt$flip)
{ results = em.shape.shift.flip(data, k=opt$cluster, iter=opt$iter, shift=opt$shift, seeding=opt$seeding)
} else { 
  results = em.shape.shift(data, k=opt$cluster, iter=opt$iter, shift=opt$shift, seeding=opt$seeding)
}

# write results
f = write.object(opt$output, results)

# exit code
# in case of failure, last element of the list is supposed to be NULL
l = length(results)
if(is.null(results[[l]]))
{ quit(status=1)
} else { 
  quit(status=0)
}


library("ChromDMM")
library("optparse")

#TRUE data
################## Enhancers 10 modifications ##########################
rm(list=ls())
opt=list()
opt$data="experiment_data/data.RDS"
opt$cluster=c(1,2,3,4,5,6,7,8)
opt$bin.size=40
opt$window=2000
opt$verbose=TRUE
opt$shift=21
opt$flip=TRUE
opt$seed.boolean=FALSE
opt$seed=FALSE
opt$eta=1.1
opt$nu=0.1
opt$etah=10
opt$nuh=10
opt$EM.max.iter=1000
opt$EM.num.tol=1e-06
opt$BFGS.max.iter=1000
opt$BFGS.num.tol=1e-06
opt$initialisation="kmeans++"
opt$soft.kmeans.maxit=1000
opt$soft.kmeans.stiffness=50
opt$soft.kmeans.randomInit=TRUE
opt$repetition=5
opt$parallel=TRUE


signals_name="1000_enhancers_bin_1_window_4000_only5prime.RData"
mod_subset1=c("CTCF", "H2AZ", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3")
mod_subset2=c("H3K79me2", "H3K9ac", "H3K9me3", "H4K20me1", "RNAPOL2", "DNase-seq", "MNase-seq")

#opt$bin.size


mod_number=10 
mod_subset=c("H2AZ", "H3K27ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K79me2", "H3K9ac", 
"RNAPOL2", "DNase-seq", "MNase-seq")

#other version
mod_number=4
mod_subset=c("H3K27ac", "H3K4me1", "RNAPOL2", "MNase-seq")
#opt$cluster1:8 

#opt$EM.max.iter=1000
N=1000

load(paste0("/Volumes/scratch/cs/csb/projects/enhancer_prediction/experiments/RProjects/preprint_devel/Data/K562/data_R/",signals_name))
signals=normalized_profiles  


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

#data given for ChromDMM
win.start <- nrow(signals[[1]])/2-opt$window/2+1
win.end <- nrow(signals[[1]])/2+opt$window/2



data.smaller.window <- function(d) {
  lapply(d, function(dat){
    newd <- matrix(0, ncol(signals[[1]]), opt$window)
    for (i in 1:ncol(signals[[1]])) {
      newd[i,] <- dat[ (win.start):(win.end),i]
    }
    newd
  })
}
signals.smaller.window=data.smaller.window(signals)

#data to be realigned
unshifting_window=opt$window+opt$bin.size*(opt$shift-1)
win.unshift.start <- nrow(signals[[1]])/2-unshifting_window/2+1
win.unshift.end <- nrow(signals[[1]])/2+unshifting_window/2

data.tobe.unshifted <- function(d) {
  lapply(d, function(dat){
    newd <- matrix(0, ncol(signals[[1]]), unshifting_window)
    for (i in 1:ncol(signals[[1]])) {
      newd[i,] <- dat[ (win.unshift.start):(win.unshift.end),i]
    }
    newd
  })
}


tobe.unshifted.data <- data.tobe.unshifted(signals)

signals.tobe.unshifted.binned <- Map( function(mod, modname) {
  count <- bin_signal(as.matrix(mod), opt$bin.size)
  #str(count)
  rownames(count) <- paste0( seq(nrow(count)))
  colnames(count) <- paste0('bin', seq(ncol(count)))
  
  count[which(count<0)]=0
  count=round(count)
  
  
  count
  
}, tobe.unshifted.data, names(tobe.unshifted.data))





save_path="experiment_data/"
figure_path="figures/"


#library("RcppGSL")
#library("ChromDMM")

# Bin signals ============================================================

signals.binned <- Map( function(mod, modname) {
  count <- bin_signal(as.matrix(mod), 1)
  #str(count)
  rownames(count) <- paste0(modname, '.', seq(nrow(count)))
  colnames(count) <- paste0('bin', seq(ncol(count)))
  
  count[which(count<0)]=0
  count=round(count)
  
  
  count
  
}, signals.smaller.window, names(signals.smaller.window))


signals.binned2 <- Map( function(mod, modname) {
  count <- bin_signal(as.matrix(mod), opt$bin.size)
  #str(count)
  rownames(count) <- paste0(modname, '.', seq(nrow(count)))
  colnames(count) <- paste0('bin', seq(ncol(count)))
  
  count[which(count<0)]=0
  count=round(count)
  
  
  count
  
}, signals.smaller.window, names(signals.smaller.window))

signals.subset <- signals.binned[mod_subset]
signals.subset2<-signals.binned2[mod_subset]



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


data=list()
data$data=signals.binned[mod_subset]
data$binned.data=signals.binned2[mod_subset]
data$tobe.unshifted.unflipped.binned.data=signals.tobe.unshifted.binned[mod_subset]
saveRDS(data, paste0("experiment_data/",N,"_enhancers_",length(mod_subset),"modifications.Rds"))




#Remove samples which are zero for some modification? almost everything get's removed?
zero.indices <- which(apply(sapply(signals.subset, rowSums), 
                            1, prod) == 0)

if (length(zero.indices) > 0) {
  warning(sprintf('Excluding are all-zero row(s) in the data: %s\n',
                  paste(names(zero.indices), collapse = ',')))
}


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


subset.fits <- ChromDMM::dmn(count=signals.subset, K=1:Kmax, 
                                                 bin.width=bin_size, 
                                                 verbose = F, eta=eta, nu=nu, etah=etah, nuh=nuh, S=S,
                                                 shift.reads = shift.reads, flip=flip, zeta=zeta, xi=xi,
                                                 repetition=repetition, randomInit=T, parallel=parallel,EM.threshold=1e-6, 
                                                 maxNumOptIter = 1000, numOptRelTol=1e-6, maxIt=maxIt, init="kmeans++" )

#histogram of the learned shift states
learned.shift.amounts=-bin.size*( floor(S/2) + 1 ) + bin.size*s_sum

dir.create(path=paste0(figure_path, gsub(".Rds", "", save_name),"_K_", cluster_nro), recursive=T)

png(paste0(figure_path, gsub(".Rds", "", save_name),"_K_", cluster_nro,
           '/learned.shift.amounts.png'), 1200, 800, res = 82)
hist(learned.shift.amounts, breaks=seq(-400,400, 40), include.lowest = TRUE, right=TRUE)
dev.off()

unshifted.binned.data=signals.tobe.unshifted.binned

#s=attr(subset.fits[[cluster_nro]]@Data, "shifts")
#learned.flip.states=attr(subset.fits[[cluster_nro]]@Data, "flips")

s=s_sum
learned.flip.states=f_cond


Lx=50
N=nrow(signals.tobe.unshifted.binned[[1]]) #can be less than 1000

#(ncol(signals.tobe.unshifted.binned[[1]]) is 70
#unshifting_window/bin.size is also 70

for(m in 1:length(unshifted.binned.data)){
  unshifted.binned.data[[m]]=matrix(0,N, ncol=Lx)
  for(i in 1:N){
    #print( seq((S-s[i]+1) ,(ncol(binned.data[[1]])-s[i]+1),1) )
    if(learned.flip.states[i]==1){
      unshifted.binned.data[[m]][i,]=signals.tobe.unshifted.binned[[m]][i, seq((S-s[i]+1) ,(ncol(signals.tobe.unshifted.binned[[1]])-s[i]+1),1)]
    }else{
      #print(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
      #length(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
      unshifted.binned.data[[m]][i,]=rev( signals.tobe.unshifted.binned[[m]][i, seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1)] )
    }
  }
  rownames(unshifted.binned.data[[m]]) <- paste0('loci', seq_len(N))
  colnames(unshifted.binned.data[[m]]) <- paste0('pos', seq_len(Lx))
}

# unshifted.binned.data2=signals.tobe.unshifted.binned

#s=attr(subset.fits[[cluster_nro]]@Data, "shifts")
#learned.flip.states=attr(subset.fits[[cluster_nro]]@Data, "flips")

# s=S-s_sum+1
# learned.flip.states=3-f_cond
# 
# 
# Lx=50
# N=nrow(signals.tobe.unshifted.binned[[1]]) #can be less than 1000
# for(m in 1:length(unshifted.binned.data2)){
#   unshifted.binned.data2[[m]]=matrix(0,N, ncol=Lx)
#   for(i in 1:N){
#     #print( seq((S-s[i]+1) ,(ncol(binned.data[[1]])-s[i]+1),1) )
#     if(learned.flip.states[i]==1){
#       unshifted.binned.data2[[m]][i,]=signals.tobe.unshifted.binned[[m]][i, seq((S-s[i]+1) ,(ncol(signals.tobe.unshifted.binned[[1]])-s[i]+1),1)]
#     }else{
#       #print(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
#       #length(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
#       unshifted.binned.data2[[m]][i,]=rev( signals.tobe.unshifted.binned[[m]][i, seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1)] )
#     }
#   }
#   rownames(unshifted.binned.data2[[m]]) <- paste0('loci', seq_len(N))
#   colnames(unshifted.binned.data2[[m]]) <- paste0('pos', seq_len(Lx))
# }


plot.heatmap.final(fits=subset.fits, show.param=TRUE, 
                   sqroot=sqroot,
                   concat=FALSE,
                   M.orig=M.orig,
                   figure_path=paste0(figure_path, gsub(".Rds", "", save_name),"_K_",cluster_nro,"/"), 
                   figure_name=paste0( "heat-better"),
                   evaluation=BIC,
                   title.margin=0.1,
                   param.height=param.height,
                   data.height=data.height,
                   param="DirichletParameters",
                   pdf_width=pdf_width,
                   pdf_height=pdf_height, #If M.orig==4 add 9
                   figure_width=2*figure_width,
                   figure_height=2*figure_height,
                   figure_res=80,
                   legend_down=legend_down,
                   legend_downshift=legend_downshift,
                   DNAlabels=DNAlabels2,
                   cluster_nro=cluster_nro, data=unshifted.binned.data, shift=TRUE, flip=TRUE)

plot.heatmap.final(fits= subset.fits, show.param=TRUE, 
                   sqroot=sqroot,
                   concat=FALSE,
                   M.orig=M.orig,
                   figure_path=paste0(figure_path, gsub(".Rds", "", save_name),"_K_",cluster_nro,"/"),
                   figure_name=paste0( "heat-better"),
                   evaluation=BIC,
                   title.margin=0.1,
                   param.height=param.height,
                   data.height=data.height,
                   param="average",
                   pdf_width=pdf_width,
                   pdf_height=pdf_height,
                   figure_width=2*figure_width,
                   figure_height=2*figure_height,
                   figure_res=80,
                   legend_down=legend_down,
                   legend_downshift=legend_downshift,
                   DNAlabels=DNAlabels2,
                   cluster_nro=cluster_nro, data=unshifted.binned.data, shift=TRUE)

print(getwd())
print(paste0(save_path, save_name))
saveRDS(subset.fits, paste0(save_path, save_name))
saveRDS(zero.indices, paste0(save_path, gsub(".Rds", "_zeroindices.Rds", save_name)))




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
              help="Print the progress of the EM algorithm"),
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






opt=list()
opt$data="experiment_data/data.RDS"
opt$cluster=c(1,2,3)
opt$bin.size=1
opt$window=50
opt$verbose=TRUE
opt$shift=21
opt$flip=TRUE
opt$seed.boolean=FALSE
opt$seed=FALSE
opt$eta=1.1
opt$nu=0.1
opt$etah=10
opt$nuh=10
opt$EM.max.iter=250
opt$EM.num.tol=1e-06
opt$BFGS.max.iter=1000
opt$BFGS.num.tol=1e-06
opt$initialisation="kmeans++"
opt$soft.kmeans.maxit=1000
opt$soft.kmeans.stiffness=50
opt$soft.kmeans.randomInit=TRUE
opt$repetition=5
opt$parallel=TRUE
opt$output="data_experiments/results.RDS"



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


library("optparse")
library("ChromDMM")
# options
option_list = list(
  make_option(c("-d", "--data"), action="store", default=NULL, type='character',
              help="The file containing the data to cluster."),
  make_option(c("--bin.size", "-B"), action="store", default=1, type='numeric',
              help="Bin size if the data resolution needs to be decreased. No change in bin size in default."),
  make_option(c("--cluster", "-c"), action="store", default=1, type='numeric',
              help="Number of clusters"),
  make_option(c("-f", "--fit"), action="store", default=NULL, type='character',
              help="The file containing the ChromDMM fit."),
  make_option(c("--fig.width"), action="store", default=1100, type='numeric',
              help="png figure width"),
  make_option(c("--fig.height"), action="store", default=1400, type='numeric',
              help="png figure height"),
  make_option(c("--name"), action="store", default="", type='character',
              help="The figure file name")
)




opt = parse_args(OptionParser(option_list=option_list))

#source("scripts/plot.heatmap.final.R")
#opt=list()
#opt$data=data
#opt$name=name
#opt$fit=fit
#opt$bin.size=1
#opt$cluster=6

#opt$data="experiment_data/data.RDS"
#opt$fit="experiment_data/simulated_data_fit_mac.RDS"
#opt$name="simulated_data_fit"

data<- readRDS( opt$data) 
window=ncol(data$data[[1]])*opt$bin.size
print(paste0("window: ", window))

fits=NULL
if( !is.null(opt$fit)){
  fits=readRDS( opt$fit)
  
  #Remove samples which are zero for some modification? almost everything get's removed?
  zero.indices <- which(apply(sapply(data$data, rowSums), 
                              1, prod) == 0)
  
  if(length(zero.indices)!=0){
    data$data=lapply(data$data, function(x)x[-zero.indices,])
    data$binned.data=lapply(data$binned.data, function(x)x[-zero.indices,])
    data$tobe.unshifted.unflipped.binned.data=lapply(data$tobe.unshifted.unflipped.binned.data, function(x)x[-zero.indices,])
  }
  
  f=attr(fits[[opt$cluster]]@Data, "flips")
  s=attr(fits[[opt$cluster]]@Data, "shifts")
  
  #maximum number of shifts
  S=dim(fits[[opt$cluster]]@Ez)[2]
  window2=ncol(data$tobe.unshifted.unflipped.binned.data[[1]])-S+1
  print(paste0("window2: ", window2))
  window=window2*opt$bin.size
  N=nrow(data$data[[1]])
  binned.data <- data$tobe.unshifted.unflipped.binned.data
  unshifted.binned.data=data$tobe.unshifted.unflipped.binned.data #to have unshifted.binned.data in correct list structure
  
  for(m in 1:length(binned.data)){
    
    unshifted.binned.data[[m]]=matrix(0, nrow=N, ncol=window2)
    for(i in 1:nrow(binned.data[[1]])){
      #    print( seq((S-s[i]+1) ,(ncol(binned.data[[1]])-s[i]+1),1) )
      if(f[i]==1){
        unshifted.binned.data[[m]][i,]=binned.data[[m]][i, seq((S-s[i]+1) ,(ncol(binned.data[[1]])-s[i]+1),1)]
      }else{
        #print(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
        #length(seq(s[i],(unshifting_window/bin.size-(S-s[i])) ,1))
        unshifted.binned.data[[m]][i,]=rev(binned.data[[m]][i, seq(s[i] ,( ncol(binned.data[[1]])-(S-s[i]) ),1)]  )
      }
      rownames(unshifted.binned.data[[m]]) <- paste0('loci', seq_len(N))
      colnames(unshifted.binned.data[[m]]) <- paste0('pos', seq_len(window2))
    }
  }
  
  plot.heatmap.final(fits=fits, show.param=TRUE, 
                                      sqroot=FALSE,
                                      concat=FALSE,
                                      figure_path="figures/", 
                                      figure_name=opt$name,
                                      evaluation=BIC,
                                      title.margin=0.1,
                                      param.height=0.4,
                                      data.height=1.3,
                                      figure_width=opt$fig.width,
                                      figure_height=opt$fig.height,
                                      figure_res=160,
                                      legend_down=0.25,
                                      legend_downshift=0.05,
                                      param="DirichletParameters",
                                      DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb")),
                                      cluster_nro=opt$cluster, data=unshifted.binned.data, shift=TRUE)
  
  #After clustering, averages
 plot.heatmap.final(fits=fits,show.param=TRUE, 
                                      sqroot=FALSE,
                                      concat=FALSE,
                                      figure_path="figures/", 
                                      figure_name=opt$name,
                                      evaluation=BIC,
                                      title.margin=0.1,
                                      param.height=0.4,
                                      data.height=1.3,
                                      param="average",
                                      figure_width=opt$fig.width,
                                      figure_height=opt$fig.height,
                                      figure_res=160,
                                      legend_down=0.25,
                                      legend_downshift=0.05,
                                      DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb")),
                                      cluster_nro=opt$cluster, data=unshifted.binned.data, shift=TRUE)
  
  
    
}else{

plot.heatmap.final(fits=fits,
                   show.param=TRUE,
                   sqroot=FALSE,
                   concat=FALSE,
                   labels=data$labels,
                   figure_path="figures/",
                   figure_name=paste0("shifted-flipped-data"),
                   evaluation=BIC,
                   title.margin=0.1,
                   param.height=0.4,
                   data.height=1.3,
                   param="average", 
                   figure_width=opt$fig.width,
                   figure_height=opt$fig.height,
                   figure_res=160,
                   legend_down=0.25,
                   legend_downshift=0.05,
                   DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb")),
                   cluster_nro=opt$cluster, data=data$data)

}

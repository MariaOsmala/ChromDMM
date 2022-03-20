library("optparse")

# options
option_list = list(
  make_option(c("-d", "--data"), action="store", default=NULL, type='character',
              help="The file containing the data to cluster."),
  make_option(c("--bin.size", "-B"), action="store", default=1, type='numeric',
              help="Bin size if the data resolution needs to be decreased. No change in bin size in default."),
  make_option(c("--bin.data"), action="store_true", default=FALSE,
              help="Whether to bin the data in data$data"),
  make_option(c("--name"), action="store", default="", type='character',
              help="The figure file name"),
  make_option(c("--fig.width"), action="store", default=700, type='numeric',
              help="png figure width"),
  make_option(c("--fig.height"), action="store", default=1000, type='numeric',
              help="png figure height")
  
  )


# opt$bin.size=40
# opt$data="experiment_data/data.RDS"
# opt$bin.data=FALSE
# opt$fig.width=700
# opt$fig.height=1000
# opt$name="shifted-flipped-data"


opt = parse_args(OptionParser(option_list=option_list))

source("scripts/plot.heatmap.before.clustering.R")
#opt=list()
#opt$data="experiment_data/data.RDS"
#opt$data="experiment_data/1000_enhancers_4modifications.Rds"
data<- readRDS( opt$data) 

if(opt$bin.data==TRUE){
  
data$data <- data$binned.data

}

window=ncol(data$data[[1]])*opt$bin.size
K=length(unique(data$labels))

plot.heatmap.before.clustering(show.param=TRUE,
                               sqroot=FALSE,
                               figure_path="figures/",
                               figure_name=paste0(opt$name),
                               evaluation=BIC,
                               title.margin=0.3,
                               param.height=0.4,
                               data.height=1.3,
                               param="average", 
                               figure_width=opt$fig.width,
                               figure_height=opt$fig.height,
                               figure_res=160,
                               legend_down=0.05,
                               legend_downshift=0.05,
                               DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb")),
                              data=data$data)




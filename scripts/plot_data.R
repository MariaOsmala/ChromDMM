library("optparse")

# options
option_list = list(
  make_option(c("-d", "--data"), action="store", default=NULL, type='character',
              help="The file containing the data to cluster."),
  make_option(c("--bin.size", "-B"), action="store", default=1, type='numeric',
              help="Bin size if the data resolution needs to be decreased. No change in bin size in default."),
 
  make_option(c("--name"), action="store", default="", type='character',
              help="The figure file name")
)





opt = parse_args(OptionParser(option_list=option_list))

source("scripts/plot.heatmap.final.R")
#opt=list()
#opt$data="experiment_data/data.RDS"
data<- readRDS( opt$data) 
window=ncol(data$data[[1]])*opt$bin.size
K=length(unique(data$labels))

plot.heatmap.final(show.param=TRUE,
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
                   figure_width=1100,
                   figure_height=1400,
                   figure_res=160,
                   legend_down=0.25,
                   legend_downshift=0.05,
                   DNAlabels=c(paste0("-",window/2000, "kb"),"0",paste0(window/2000, "kb")),
                   cluster_nro=K, data=data$data)



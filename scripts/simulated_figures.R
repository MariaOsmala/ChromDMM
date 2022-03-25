library("optparse")
library("ChromDMM")
library("ggplot2")
# options
option_list = list(
  make_option(c("-f", "--fit"), action="store", default=NULL, type='character',
              help="The results file."),
  make_option(c("--skip"), action="store", default=1, type='numeric',
              help="How many iterations to skip from the beginning"),
  make_option(c("-c", "--cluster"), action="store", default=NULL, type='character',
              help="The number of clusters."),
  make_option(c("--name"), action="store", default="", type='character',
              help="The figure file name")
)





opt = parse_args(OptionParser(option_list=option_list))


#opt=list()
#opt$fit="experiment_data/simulated_data_fit_rep1.RDS"
#fit=readRDS("experiment_data/results.RDS") #length 1
#fit=readRDS("experiment_data/simulated_data_fit_mac.RDS") #list of length 3
#fit=readRDS("experiment_data/simulated_data_fit_rep1.RDS") #list of length 3
#fit=readRDS("experiment_data/simulated_data_fit_sbatch_K_2.RDS") #length 1

figure_path="figures/"

fit=readRDS(opt$fit)

png(paste0(figure_path,"EM.diagnostics-",opt$name,".png"), 2000, 4000, res = 150)
ChromDMM::plot.EM(fit=fit[[opt$cluster]], smoothness.scale='free', skip=0, plot.det=FALSE)
dev.off()

#ChromDMM::plot.EM(fit=fit, smoothness.scale='free', skip=0, plot.det=FALSE)

#
png(paste0(figure_path,"NLL-",opt$name,".png"), 800, 800, res = 150)
gg=ggplot(fit[[opt$cluster]]@nll.data[opt$skip:nrow(fit[[opt$cluster]]@nll.data),],aes(iter-1, nll)) +
  geom_line() +geom_point(shape=1)+xlab("EM iter")+ylab("Negative log posterior")
plot(gg)
dev.off()


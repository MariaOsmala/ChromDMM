library("optparse")
library("ChromDMM")
# options
option_list = list(
  make_option(c("-f", "--fit"), action="store", default=NULL, type='character',
              help="The results file."),
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


fit=readRDS(opt$fit)
print(length(fit))

figure_path="figures/"

png(paste0(figure_path, 'BIC-',opt$name,'.png'), 1200, 1200, res = 150)
plot(sapply(fit, BIC), type='b', xlab='Number of clusters', ylab='BIC')
dev.off()

png(paste0(figure_path, 'AIC-',opt$name,'.png'), 1200, 1200, res = 150)
plot(sapply(fit, AIC), type='b', xlab='Number of clusters', ylab='AIC')
dev.off()
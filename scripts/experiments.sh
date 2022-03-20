data="experiment_data/data.RDS"
#experiment_data/data.RDS is an example of shifted and flipped simulated data that contains two chromatin features H3K4me1 and RNAPOL2 and two clusters
#The R object in data.RDS is a list of 6 with the following elements
# $ data                         :List of 2
#  ..$ H3K4me1: num [1:1000, 1:50] 3 0 1 0 0 1 2 0 0 0 ...
#  ..$ RNAPOL2: num [1:1000, 1:50] 2 0 0 0 0 0 1 1 0 0 ...
# $ labels                       : int [1:1000] 1 1 1 1 1 1 1 1 1 1 ... 
# $ true.flips                   : num [1:1000] 1 2 1 1 2 1 1 1 2 2 ...
# $ true.shift.ind               : num [1:1000] 11 13 11 10 16 15 8 12 14 6 ...
# $ true.shift.bases             : num [1:1000] 0 80 0 -40 200 160 -120 40 120 -200 ...
# $ tobe.unshifted.unflipped.data:List of 2
#  ..$ H3K4me1: num [1:1000, 1:70] 0 0 0 1 0 1 0 2 0 0 ...
#  ..$ RNAPOL2: num [1:1000, 1:70] 0 0 0 1 0 1 0 0 0 0 ...

#There are N=1000 genomic regions in this data. 
#$data contain the chromatin feature coverage counts. The original genomic window is W=2000 a
#Labels contain the true clusters, 
#the first 500 regions belong to cluster 1 and the rest belong to cluster 2
#true.flips contain the simulated flips
#true.shift.ind contain the simulated shift states, the number of shift states is 21, and the original bin.size is 40 
#-> the maximum shift is (21-1)/2*40 -400bp +400bp
#true.shift.bases correspond to simulated shifts in nucleotides
#$tobe.unshifted.unflipped.data contains the same shifted and flipped data but in a larger genomic window W=2800. 
#This data is realigned after the model inference

bin_size=1
window=50

Rscript scripts/run_ChromDMM.R  --data $data --cluster 1,2,3 --bin.size $bin_size --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 1 --parallel TRUE --output "experiment_data/simulated_data_fit_mac.RDS"


#true enhancer data

#To see how the 1000_enhancers_10modifications.Rds and 1000_enhancers_4modifications.Rds were created from experiments_data/.Rds 
#which was in turn created from ENCODE K562 data with PREPRINT preprocessing steps by setting 
#  five_prime_end: TRUE in config.yaml file

data="experiment_data/1000_enhancers_10modifications.Rds"
#A list of 3
#$ data is a list of 10. Each list element contains N x W chromatin feature data matrices. N is the number of elements and W is the genomic window
#This data is given to the clustering
#$ binned.data, same as data, but binned using bin.size=40
#$ tobe.unshifted.unflipped.binned.data:List of 10
#Data extracted originally in window 2800 and binned 

bin_size=40
window=2000
7*5

RScript scripts/run_ChromDMM.R  --data $data --cluster 2,3,4,5,6,7,8 --bin.size $bin_size --window $window --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 5 --parallel TRUE --output "data_experiments/10mods_fit.RDS"

data="experiment_data/1000_enhancers_4modifications.Rds"
bin_size=40
window=2000

RScript scripts/run_ChromDMM.R  --data $data --cluster 6 --bin.size $bin_size --window $window --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 4 --parallel TRUE --output "data_experiments/4mods_fit.RDS"


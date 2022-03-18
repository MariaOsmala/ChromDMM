#!/bin/bash
#SBATCH --time=2-00:00:00    #2 day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    #8
#SBATCH --mem-per-cpu=5G     #10G		
#SBATCH -J 10mods
#SBATCH -o 10mods.out 
#SBATCH -e 10mods.err







module load anaconda

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/

source activate /m/cs/scratch/csb/projects/enhancer_clustering/EnhancerClustering/experiments/MariasExperimentsPackageProperlyInstalled/enhancerdnase-experiments/DMM_experiments


export OMP_NUM_THREADS=24

source activate ./ChromDMM


srun Rscript --vanilla run_dmn_only5prime.R 


data="experiment_data/1000_enhancers_10modifications.Rds"


bin_size=40
window=2000

srun Rscript scripts/run_ChromDMM.R  --data $data --cluster 2,3,4,5,6,7,8 --bin.size $bin_size --window $window --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 5 --parallel TRUE --output "data_experiments/10mods_fit.RDS"



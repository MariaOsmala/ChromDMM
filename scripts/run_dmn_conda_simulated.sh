#!/bin/bash
#SBATCH --time=12:00:00    #2 day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12    #8
#SBATCH --mem-per-cpu=5G     #10G		
#SBATCH -J simulated
#SBATCH -o simulated.out 
#SBATCH -e simulated.err







module load anaconda

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/

source activate /m/cs/scratch/csb/projects/enhancer_clustering/EnhancerClustering/experiments/MariasExperimentsPackageProperlyInstalled/enhancerdnase-experiments/DMM_experiments


export OMP_NUM_THREADS=12

source activate ./ChromDMM


srun Rscript --vanilla run_dmn_only5prime.R 



data="experiment_data/data.RDS"

bin_size=1
window=50

srun Rscript scripts/run_ChromDMM.R  --data $data --cluster 1,2,3 --bin.size $bin_size --window $window --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 5 --parallel TRUE --output "data_experiments/simulated_data_fit_sbatch.RDS"


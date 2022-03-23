#!/bin/bash
#SBATCH --time=2-00:00:00    #2 day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12    #8
#SBATCH --mem-per-cpu=1G     #10G		
#SBATCH -J simulated
#SBATCH -o simulated.out 
#SBATCH -e simulated.err







module load anaconda

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/


export OMP_NUM_THREADS=12

source activate ./ChromDMM


data="experiment_data/data.RDS"

bin_size=1
window=50

#This worked with repetition 1

srun Rscript scripts/run_ChromDMM.R  --data $data --cluster 1,2,3 --bin.size $bin_size --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 10 --parallel TRUE --output "experiment_data/simulated_data_fit_batch_10reps.RDS"


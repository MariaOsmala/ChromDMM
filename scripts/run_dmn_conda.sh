#!/bin/bash
#SBATCH --time=5-00:00:00    #2 day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24   #8
#SBATCH --mem-per-cpu=3G     #10G		
#SBATCH -J 10mods_10reps
#SBATCH -o 10mods_10reps.out 
#SBATCH -e 10mods_10reps.err







module load anaconda

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/



export OMP_NUM_THREADS=24

source activate ./ChromDMM




data="experiment_data/1000_enhancers_10modifications.Rds"


bin_size=40
window=2000

srun Rscript scripts/run_ChromDMM.R  --data $data --cluster 2,3,4,5,6,7,8 --bin.size $bin_size --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 10 --parallel TRUE --output "experiment_data/10mods_fit_10_rep.RDS"



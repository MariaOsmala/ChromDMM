#!/bin/bash
#SBATCH --time=2-00:00:00    #2 day #with 1 replication this took 9 hours
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24    #8
#SBATCH --mem-per-cpu=1G     #10G		
#SBATCH -J 4mods
#SBATCH -o 4mods.out 
#SBATCH -e 4mods.err







module load anaconda

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/


export OMP_NUM_THREADS=24

source activate ./ChromDMM




data="experiment_data/1000_enhancers_4modifications.Rds"


bin_size=40
window=2000

srun Rscript scripts/run_ChromDMM.R  --data $data --cluster 2,3,4,5,6,7,8 --bin.size $bin_size --verbose FALSE --shift 21 --flip TRUE --seed.boolean FALSE --repetition 10 --parallel TRUE --output "experiment_data/4mods_fit_10_rep.RDS"



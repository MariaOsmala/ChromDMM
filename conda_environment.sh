These packages are required to install ChromDMM, ChromDMM installs the required packages itself
- r-base
- r-essentials
- r-devtools

These are required in a linux environemnt, gcc version is crucial
  - gcc=11.2
  - gxx
  - gfortran


conda env create -f conda_environment.yml
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

conda env update --file conda_environment.yml



########### In triton ###################

sinteractive --time=04:00:00 --mem-per-cpu=5G --nodes=1 --ntasks=1 --cpus-per-task=15

cd /m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/
module load anaconda


export OMP_NUM_THREADS=15

#MAMBA WOULD BE FASTER
#source activate /m/cs/scratch/csb/projects/enhancer_clustering/snakemake/mamba 

#mamba env create --prefix ./ChromDMM_mamba --force -f conda_environment.yml
#mamba config --set env_prompt '({name})'
#source activate ./preprint_mamba
#conda config --env --add channels defaults
#conda config --env --add channels bioconda
#conda config --env --add channels conda-forge


conda env create --prefix ./ChromDMM --force -f conda_environment.yml
conda config --set env_prompt '({name})'
source activate ./ChromDMM

Rscript -e 'install.packages("/scratch/cs/csb/projects/enhancer_clustering/Rpackages/ChromDMM_1.0.tar.gz", repos = NULL, type = .Platform$pkgType, lib="/m/cs/scratch/csb/projects/enhancer_clustering/Rpackages/ChromDMM/ChromDMM/lib/R/library")'

conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

R
Rcpp::compileAttributes()
devtools::install()





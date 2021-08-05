
#update conda by running
Conda update -n base -c defaults conda

#macillä
Luo kotihakemistoon .zshrc tiedosto ja kirjoita siihen
export PATH=$HOME/anaconda3/bin:$PATH

#Then run
vi ~/.zshrc 
. ~/.zshrc 

conda #it works!

conda info --envs
conda env list


conda env create -f DirichletMultinomial_R4.0.3.yml

#Overrride
conda env create --force -f experiments/DMM_experiments_R4.0.3.yml 

source activate DirichletMultinomial_R4.0.3

conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

#Tämän jälkeen
gsl-config --cflags #-I/Users/mpirttin/anaconda3/envs/DirichletMultinomial/include
gsl-config --libs #-L/Users/mpirttin/anaconda3/envs/DirichletMultinomial/lib -lgsl -lgslcblas

export CFLAGS="-I/Users/mpirttin/anaconda3/envs/DirichletMultinomial_R4.0.3/include"
export LDFLAGS="-L/Users/mpirttin/anaconda3/envs/DirichletMultinomial_R4.0.3/lib -lgsl -lgslcblas"

export RSTUDIO_WHICH_R=/Users/mpirttin/anaconda3/envs/DirichletMultinomial_R4.0.3/bin/R

	open -na Rstudio




#In RStudio
.libPaths("/Users/mpirttin/anaconda3/envs/DirichletMultinomial_R4.0.3/lib/R/library")

#Update package

The --prune option causes conda to remove any dependencies that are no longer required from the environment.

conda env update DirichletMultinomial_R4.0.3 --file DirichletMultinomial_R4.0.3.yml  --prune


make an exact copy of an environment by creating a clone of it:
conda create --name myclone --clone myenv

conda list --explicit

To use the spec file to create an identical environment on the same machine or another machine:

conda create --name myenv --file spec-file.txt
To use the spec file to install its listed packages into an existing environment:

conda install --name myenv --file spec-file.txt 

Export your active environment to a new file:

conda env export --from-history > environment.yml

Create the environment from the environment.yml file:

conda env create -f environment.yml




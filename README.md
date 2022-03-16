# ChromDMM

This R package contains ChromDMM, a Dirichlet-multinomial mixture model for clustering heterogeneous epigenetic data

## Installation

You can install the R package using [devtools](https://devtools.r-lib.org):

```R
devtools::install_github('MariaOsmala/ChromDMM')
```

ChromDMM depends on the following R-packages, which should be installed automatically when you install the `ChromDMM` package: CHECK THESE!

  - [bioconductor](http://bioconductor.org)
  - [bioconductor-rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
  - [bioconductor-genomeinfodb](http://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html)
  - [bioconductor-genomeinfodbdata](http://bioconductor.org/packages/release/data/annotation/html/GenomeInfoDbData.html)
  - [bioconductor-bsgenome.hsapiens.ucsc.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)
  - [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
  - [doparallel](https://cran.r-project.org/web/packages/doParallel/index.html)
  - [gdata](https://cran.r-project.org/web/packages/gdata/index.html)
  - [MASS](https://cran.r-project.org/web/packages/MASS/index.html)
  - [caret](https://cran.r-project.org/web/packages/caret/index.html)
  - [e1071](https://cran.r-project.org/web/packages/e1071/index.html)


## Running (part of) the analysis pipeline in the paper

The scripts to reproduce a part of the analysis pipeline as described in the paper, are in the `ChromDMM/scripts` folder.

Running the analysis depends on more software than just the ChromDMM package.
To install the required software, R and Python packages, we recommend using [anaconda](https://www.anaconda.com/products/individual) and the provided `conda_environment.yml` file:

```
conda env create -f ChromDMM/conda_environment.yml -n ChromDMM
source activate ChromDMM
```

Alternatively, you could install the software manually.

Software:

  - [R](https://www.r-project.org) (version 3)
  - [Python](https://python.org)
  - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - [samtools](http://www.htslib.org)
  - [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
  - [phantompeakqualtools](https://www.encodeproject.org/software/phantompeakqualtools/)

Python packages:

  - [snakemake](https://snakemake.readthedocs.io/en/stable/)
  - [pandas](https://pandas.pydata.org)

When all required software has been installed, you can continue to run the analysis by first editing `preprint/workflows/config.yaml` to specify the folders in which to place the data (needs around 1 terabyte of free space).
Then, from the main `preprint/` folder, run:

```
snakemake
```

This downloads all data and runs the entire analysis pipeline.

## Citation:

Osmala, M., Eraslan, G. & Lähdesmäki, H. (2022). ChromDMM: A Dirichlet-Multinomial Mixture Model For Clustering Heterogeneous Epigenetic Data (submitted.)
## License

This project is licensed under the LGPL-3 License - see the [LICENSE.md](LICENSE.md) file for details

## Contact

Maria Osmala, MSc  
PhD student  
Aalto University School of Science  
Department of Computer Science  
Email: firstname.surname@aalto.fi  
Home Page: https://people.aalto.fi/maria.osmala

Harri Lähdesmäki, D. Sc. (Tech)  
Associate Professor  
Aalto University School of Science  
Department of Computer Science  
Email: firstname.surname@aalto.fi  
Home Page: http://users.ics.aalto.fi/harrila

## Acknowledgments

* 
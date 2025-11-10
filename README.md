## To get started using conda
```{bash}
conda config --add channels conda-forge
conda config --add channels bioconda

conda create --prefix ${YOUR_ENV_PATH} r-base=4.4 r-devtools r-seurat -c conda-forge
conda activate ${YOUR_ENV_PATH}
conda install r-signac -c bioconda
```
- You'll need R version >= 4.4 to run the package.

## Install geneSetSimplifyR and fGSEA in R
```{r}
library(devtools)
install_github("alserglab/fgsea")
install_github("chihs-dtu/geneSetSimplifyR", build_vignettes=TRUE)
```

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
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fgsea")

devtools::install_github("chihs-dtu/geneSetSimplifyR", build_vignettes=TRUE)
```
- If you cannot install the package through `install_github`, clone the respository manually then install the package.
```{bash}
# Run this in bash
git clone git@github.com:chihs-dtu/geneSetSimplifyR.git

# Run this in R with the directory where the repository is cloned.
devtools::install("PATH_TO_DIRECTORY")
```

## Start using geneSetSimplifyR
- check the vignette [here](vignettes/geneSetSimplifyR.html)


## Install geneSetSimplifyR and fGSEA in R
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fgsea")

devtools::install_github("chihs-dtu/geneSetSimplifyR", build_vignettes=TRUE)
```
If you are using Apple Silicon and encounter the error message "ERROR: sub-architecture 'R' is not installed", add `INSTALL_opts = c("--no-multiarch")` to both `BiocManager::install()` and `devtools::install_github()`.

## Start using geneSetSimplifyR
- check the vignette [here](vignettes/geneSetSimplifyR.html)


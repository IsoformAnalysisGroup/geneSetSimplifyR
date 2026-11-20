## Creates data/exampleGeneSets.rda, the gene sets of the example.
##
## 'exampleGeneSets' holds the human MSigDB v7.5.1 gene sets (gene symbols)
## of 'exampleEnrichment', from the hallmark (H), curated (C2: CGP, Reactome,
## WikiPathways) and ontology (C5: GO) collections. Each gene set is
## restricted to the genes tested in the TCGA colon adenocarcinoma (COAD)
## analysis of Vitting-Seerup & Sandelin (2017),
## https://doi.org/10.1158/1541-7786.MCR-16-0459. KEGG gene sets are left out
## because the KEGG license does not allow their redistribution.
##
## The list of tested genes was not kept. The genes of the shipped gene sets
## are exactly the tested genes in those gene sets, so restricting the MSigDB
## gene sets to them gives the same result.
##
## Needs msigdbr 7.5.1 (CRAN archive), the last version with MSigDB v7.5.1.
## Run from the root of the package source:
##   Rscript inst/scripts/make-exampleGeneSets.R
## and then inst/scripts/make-exampleEnrichment.R and make-exampleGsList.R.

stopifnot(packageVersion("msigdbr") == "7.5.1")

load(file.path("data", "exampleGeneSets.rda"))
testedGenes <- unique(unlist(exampleGeneSets))

msigdb <- msigdbr::msigdbr(species = "Homo sapiens")
msigdb <- msigdb[msigdb$gs_name %in% names(exampleGeneSets) &
    msigdb$gs_subcat != "CP:KEGG", ]
keep <- intersect(names(exampleGeneSets), msigdb$gs_name)

exampleGeneSets <- lapply(
    split(msigdb$gene_symbol, msigdb$gs_name)[keep],
    function(genes) intersect(genes, testedGenes)
)

save(
    exampleGeneSets,
    file = file.path("data", "exampleGeneSets.rda"),
    compress = "xz"
)

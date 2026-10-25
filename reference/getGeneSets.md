# Fetch MSigDB gene sets by pathway name

Fetches gene-set memberships from MSigDB via
[`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
and returns them in the named-list format expected by
[`geneSetSimplifyR()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/geneSetSimplifyR.md)
(pathway name -\> character vector of gene IDs).

Useful in two situations:

- Inspecting the gene members of specific pathways (e.g. the top
  pathways from a cluster after running the workflow).

- Building a `geneSetsList` from scratch without going through
  `pairedGSEA::prepare_msigdb()`.

## Usage

``` r
getGeneSets(
  pathway_names = NULL,
  species = "Homo sapiens",
  collection = "C5",
  subcollection = NULL,
  gene_id_type = c("ensembl_gene", "gene_symbol", "entrez_gene")
)
```

## Arguments

- pathway_names:

  Character vector of MSigDB pathway names (i.e. `gs_name` values). If
  `NULL` (default), all pathways in the requested
  collection/subcollection are returned.

- species:

  Species, passed to
  [`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html).
  Default `"Homo sapiens"`.

- collection:

  MSigDB collection (e.g. `"C5"`, `"H"`, `"C2"`). Default `"C5"`.

- subcollection:

  Optional MSigDB sub-collection (e.g. `"GO:BP"`). `NULL` (default)
  returns the full collection.

- gene_id_type:

  Which gene identifier to use for the returned vectors. One of
  `"ensembl_gene"` (default), `"gene_symbol"`, or `"entrez_gene"`.

## Value

A named list where each element is a character vector of gene IDs for a
pathway. Names are pathway names (`gs_name`).

## Examples

``` r
# \donttest{
if (requireNamespace("msigdbr", quietly = TRUE)) {
  # Hallmark gene sets (gene symbols)
  hallmark <- getGeneSets(collection = "H", gene_id_type = "gene_symbol")

  # Just the genes of two specific pathways
  twoSets <- getGeneSets(
    pathway_names = c("HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA"),
    collection = "H",
    gene_id_type = "gene_symbol"
  )
  lengths(twoSets)
}
#> Downloading gene sets (first use only, may take a few minutes)...
#> HALLMARK_APOPTOSIS   HALLMARK_HYPOXIA 
#>                161                200 
# }
```

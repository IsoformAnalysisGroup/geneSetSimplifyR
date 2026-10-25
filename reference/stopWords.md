# Stop Words

Words ignored when labeling clusters, e.g. "of", "and", "process".
Bigrams (two-word label candidates) containing one of these words are
not used as labels.

## Usage

``` r
stopWords
```

## Format

A character vector.

## Examples

``` r
head(stopWords)
#> [1] "gobp" "gomf" "hp"   "gocc" "kegg" "gcnp"
# Add your own stop words
myStopWords <- c(stopWords, "regulation", "response")
```

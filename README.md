## IMPACD Dependencies
**IMPACD** relies on a few packages from Bioconductor and github. Here is how to install them from within RStudio or the R console.
``` r
# Install bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowWorkspace", "CytoML", "flowCore"))

# Install github dependencies
devtools::install_github("RGLab/cytoinstaller")
cytoinstaller::install_cyto(bioc_ver = "devel")
devtools::install_github("DillonHammill/CytoExploreR")
```

After this, you can install **IMPACD** from github using the following command:
``` r
devtools::install_github("BenaroyaResearch/IMPACD")
```

## Usage
**IMPACD** ships with a vignette that walks through the entire analysis pipeline. You can access it by running the following command:
```r
vignette("IMPACD")
```

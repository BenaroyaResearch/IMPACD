# Installing IMPACD
## Dependencies
**IMPACD** relies on a few packages from Bioconductor and github. Here is how to install them from within RStudio or the R console.
``` r
# Install bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowWorkspace", "CytoML", "flowCore"))

# Install github dependencies
devtools::install_github("DillonHammill/CytoExploreR")
```
## Access token

The first step of installation is to generate an access token so that R knows you have permission to install this package. To do so, go to https://github.com/settings/tokens and click "Generate new token". Name the token something that makes sense to you (like "BRI repo access token") and click the "Repo" checkbox to define the token scope. It's important to keep this token private. This scope will provide full access to the repo sourcing this package, as well as any other private repos you have access to.

Copy this token somewhere you can refer back to, as you will need to use this token for the next step.

## Install from github
Once you have your token open an R session and execute the following to load your token:
```r
Sys.setenv(
 GITHUB_PAT = "OUR_TOKEN_GOES_HERE",
 R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true"
)
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

## miMediation

[miMediation](https://github.com/KiRinHong/miMediation) is an package for performing phylogeny-based mediation test [PhyloMed](https://github.com/KiRinHong/PhyloMed) for microbiome data.

See the following:

- Hong, Q., Chen G., and Tang Z-Z.. (2021) A phylogeny-based test of mediation effect in microbiome. Manuscript.

## Author

Qilin Hong @[Tang lab](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## Installation

You can install the package from github with:

``` r
devtools::install_github("KiRinHong/miMediation")
```
The Package source can be downloaded from https://tangzheng1.github.io/tanglab/software, which can be installed with:

``` r
install.packages(miMediation_0.1.tar.gz, repos = NULL, type ="source") 
```

## Vignette

If you install the package from source, a vignette describing the use of the package is available from within R. Load the package and then use the vignette function.

``` r
library(miMediation)
vignette("miMediation",package="miMediation")
```
Otherwise, it will not build vignettes by default if you install the package from github because they’re time consuming and may require additional packages. You can force building with:

``` r
devtools::install_github("KiRinHong/miMediation", build_vignettes = TRUE)
```
Then, the vignette would be availabe from within R.

## Getting help

Please use the issues tab (https://github.com/KiRinHong/miMediation/issues) to file any bugs or suggestions.

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3


<img src="inst/webapp/www/a2tea_hexsticker.png" align="right" alt="" width="120" />

<!-- README.md is generated from README.Rmd. Please edit that file -->

# A2TEA.WebApp

<!-- badges: start -->

[![R-CMD-check](https://github.com/tgstoecker/A2TEA.WebApp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tgstoecker/A2TEA.WebApp/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/github/last-commit/tgstoecker/A2TEA.WebApp.svg)](https://github.com/tgstoecker/A2TEA.WebApp/commits/master)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

A2TEA is a software workflow facilitating identification of candidate
genes for stress adaptation based on comparative genomics and
transcriptomics. It combines differential gene expression with gene
family expansion as an indicator for the evolution of adaptive traits.  
The goal of the A2TEA.WebApp is to allow exploration, highlighting, and
exporting of the results generated by a A2TEA pipeline run. The core of
the package is a Shiny application that aims to combine phylogenetic,
transcriptomic & functional data and results, in a way that makes it
easier to generate insightful observations. As the A2TEA pipeline
requires the user to formulate hypotheses - for each of which specific
results are produced - the WebApp allows to investigate the genomic
adaptations of one or a subset of the investigated species to the
analyzed treatment condition. The A2TEA.WebApp aims to combine the
benefits of interactivity and reproducibility, e.g. by capturing genes
or orthologous groups of interest bookmarked by the user and allowing
for easy export of plots, tables and complete subsets of the original
input data.

<br>

## :arrow_double_down: Download & Setup

There are currently 3 ways of using the A2TEA.WebApp:  
- installing the development version as an R package from this
repository with devtools/remotes  
- a docker container of the latest stable release, for which we offer a
singularity based guide  
- a demo instance of the App running hosted on
[shinyapps.io](https://www.shinyapps.io/) allowing for a sneak peak of
the functionality

### Option 1: Installation

You can install the development version of A2TEA.WebApp from GitHub with
either remotes or devtools, e.g.:

``` r
library(remotes)
remotes::install_github("tgstoecker/A2TEA.WebApp", 
                        dependencies = TRUE, 
                        build_vignettes = TRUE)
```

Note that some system dependencies might have to be installed:

``` bash
apt-get update
apt-get upgrade
apt-get install -y \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libcurl4-openssl-dev \
    unixodbc-dev \
    libproj-dev xdg-utils \
    --fix-missing
```

### Option 2 - Docker container using singularity

You can also circumvent dependency issues by downloading the latest
release of our [A2TEA Docker
image](https://hub.docker.com/repository/docker/tgstoecker/a2tea_webapp).  
To circumvent most potential problems regarding access rights we use
Singularity to pull and use the image.  
Note however, that we require `--writable` to be enabled and therefore
most likely require either the `--fakeroot` or `--no-home` flag as well.
Which command works will depend on your Singularity installation -
whether or not it was priviliged or not.  
If you have sudo rights and installed Singularity via your systems
package manager then there should be no problem at all.

``` bash
#pull the image from dockerhub
singularity pull a2tea_webapp.sif docker://tgstoecker/a2tea_webapp:latest

#open R console of image in non-persistent but writable mode
singularity run --writable --fakeroot a2tea_webapp.sif R
#or
singularity run --writable --no-home a2tea_webapp.sif R
```

#### Starting the App

For both option 1 and 2 starting the Shiny Application requires loading
the library and calling the A2TEA_App function:

``` r
library(A2TEA.WebApp)
A2TEA_App()
#upload an A2TEA.RData file
```

### Option 3 - Hosted version on shinyapps.io

Click on this link - <https://tgstoecker.shinyapps.io/A2TEA-WebApp/> -
and a new browser tab will open in which a running instance of the
WebApp will be started.

<br>

## Usage overview

All outlined setup options of the App come with an included demo dataset
which can be loaded via a single button press. You can find the rendered
version of the documentation of `A2TEA.WebApp` at the project website
<https://tgstoecker.github.io/A2TEA.WebApp>, created with `pkgdown`,
which includes a usage guide detailing all of the App’s functionality.

<br>

## Development

If you encounter a bug, have a usage questions, or want to share ideas
and functionality to make this package better, feel free to file an
[issue](https://github.com/tgstoecker/A2TEA.WebApp/issues).

<br>

## License

MIT © Tyll Stöcker

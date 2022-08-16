#' @import shiny
#' @import dplyr
#' @import ggtree
#' @import ape
#' @import readr
#' @import stringr
#' @import tidyr
#' @import purrr
#' @import treeio
#' @import ggtreeExtra
#' @import ggstar
#' @import ggnewscale
#' @import ggplot2
#' @import ggstance
#' @import Biostrings
#' @import ggmsa
#' @import seqmagick
#' @import XVector
#' @import gtable
#' @import grid
#' @import DT
#' @import cowplot
#' @import ggplotify
#' @import shinydashboard
#' @import shinydashboardPlus
#' @import fontawesome
#' @import UpSetR
#' @import ggvenn
#' @import shinyjs
#' @import Polychrome
#' @import shinycssloaders
#' @import svglite
#' @import topGO
#' @import Rgraphviz
#' @import scales
#' @import methods


.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("A2TEA.WebApp", 
                         system.file("www", 
                                     package = "A2TEA.WebApp"))
}


#####
# defining the A2TEA HYPOTHESES object
#####
# class for the expanded_OG - 
#containing all different types of data we have on it
setClass("expanded_OG", slots=list(blast_table="tbl_df",
                                   add_OG_analysis="list"))

# class for the hypotheses
setClass("hypothesis", slots=list(description="character", 
                                  number="character",
                                  expanded_in ="character", 
                                  compared_to="character", 
                                  expanded_OGs="list",
                                  species_tree="phylo"))

#another class for adding OGs analysis
setClass("add_OG_set",
         slots=list(genes="spec_tbl_df",
                    msa="AAStringSet", 
                    tree="phylo"
         )
)


#####
#Helper functions
#####
#hyperlinks
createNCBILink <- function(val) {
  sprintf('<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=%s" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
}

createEnsemblPlantsLink <- function(val, ensembl_db) {
  if (ensembl_db == "Plants") {
    sprintf('<a href="https://plants.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  } else if (ensembl_db == "Bacteria") {
    sprintf('<a href="https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  } else if (ensembl_db == "Protists") {
    sprintf('<a href="https://protists.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  } else if (ensembl_db == "Metazoa") {
    sprintf('<a href="https://metazoa.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  } else if (ensembl_db == "Fungi") {
    sprintf('<a href="https://fungi.ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  } else if (ensembl_db == "Vertebrates") {
    sprintf('<a href="https://ensembl.org/Multi/Search/Results?species=all;idx=;q=%s;site=ensemblunit" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
  }
}

#AmiGO link
createAmiGOLink <- function(go_term) {
  sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary btn-custom-link">%s</a>', go_term, go_term)
}

#since we might run into a problem if two user login at the very exact time;
#we can add an additional random string to the end of the session directory
#https://stackoverflow.com/questions/42734547/generating-random-strings
random_token_generator <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#FUNCTION - new directory
fn_dir <- function(currwd=NULL)
{
  if(is.null(currwd)) stop("Argument 'currwd' is empty.\n")
  
  #Create new working directory
  newwd <- paste0(currwd, "/",
                  format(Sys.time(), "%Y%m%d%H%M%S"), "_",
                  random_token_generator(1))
  dir.create(newwd)
  newwd
}
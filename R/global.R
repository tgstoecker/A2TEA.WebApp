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
                         system.file("webapp/www", 
                                     package = "A2TEA.WebApp"))
}


#####
# defining the A2TEA HYPOTHESES object
#####

#' expanded_OG Class
#' 
#' A class for expanded OG's
#' 
#' @name expanded_OG
#' @slot blast_table df. gene level; blast hits of OG and all closest OGs
#' @slot add_OG_analysis list. lists subclass add_OG_set elements
#' 
#' @export expanded_OG
expanded_OG <- setClass("expanded_OG", slots=list(blast_table="tbl_df",
                                   add_OG_analysis="list"))

#' hypothesis Class
#' 
#' A class for the hypotheses
#' 
#' @name hypothesis
#' @slot description character. name given to hypothesis
#' @slot number integer. index number of hypothesis
#' @slot expanded_in character. species investigated for expansion
#' @slot compared_to character. species that are used for comparison
#' @slot expanded_OGs list. lists all expanded OGs for current hypothesis
#' @slot species_tree phylo. species tree for current hypothesis
#' 
#' @export hypothesis
hypothesis <- setClass("hypothesis", slots=list(description="character", 
                                  number="character",
                                  expanded_in ="character", 
                                  compared_to="character", 
                                  expanded_OGs="list",
                                  species_tree="phylo"))

#' add_OG_set Class
#' 
#' A class for the subelements of the add-OGs analysis
#' 
#' @name add_OG_set
#' @slot genes spec_tbl_df. DF of genes for current set
#' @slot msa AAStringSet. MSA for current set
#' @slot tree phylo. Phylogenetic tree of current set
#' 
#' @export add_OG_set
add_OG_set <- setClass("add_OG_set",
         slots=list(genes="spec_tbl_df",
                    msa="AAStringSet", 
                    tree="phylo"
         )
)


#####
#Helper functions
#####

#credit to the GeneTonic devs for this nice add-on for DT log2FC display
#' Adding barplot overlay to DataTables
#'
#' @param data Column choice of underlying dataframe which shoud receive barplot overlay.
#' @param color_pos Color choice for positive values.
#' @param color_neg Color choice for negative values.
#'
#' @return Add styling based on choices.
#'
#' @export
styleColorBar_divergent <- function(data,
                                    color_pos,
                                    color_neg) {
  
  max_val <- max(abs(data))
  JS(
    sprintf(
      "isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
      max_val, color_pos, max_val, color_pos, color_neg, color_neg, max_val, max_val))
}

#AmiGO link
#' Create NCBI link
#' 
#' @param val Gene of interest
#' 
#' @return Return NCBI link for gene of interest
#' 
#' @export
createNCBILink <- function(val) {
  sprintf('<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=%s" target="_blank" class="btn btn-primary btn-custom-link">Link</a>',val)
}

#' Create Ensembl link
#'
#' @param val Gene of interest
#' @param ensembl_db Which ensembl database to use
#' 
#' @return Return NCBI link for gene of interest
#' 
#' @export
createEnsemblLink <- function(val, ensembl_db) {
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
#' Create AmiGO link
#' 
#' @param go_term Go term of interest
#' 
#' @return Return AmiGO link for go_term
#' 
#' @export
createAmiGOLink <- function(go_term) {
  sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary btn-custom-link">%s</a>', go_term, go_term)
}

#since we might run into a problem if two user login at the very exact time;
#we can add an additional random string to the end of the session directory
#https://stackoverflow.com/questions/42734547/generating-random-strings

#' Create a random token
#' 
#' @param n sample size
#' @return New random token is created
#' 
#' @export
random_token_generator <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#' Create a new tmp directory for plots, etc.
#' 
#' @param currwd Current working directory.
#' @return Creates a new tmp directory that is deleted after exiting session.
#' 
#' @export
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
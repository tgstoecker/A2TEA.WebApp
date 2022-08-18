#!/usr/bin/env r

remotes::install_github('tgstoecker/A2TEA.WebApp'); 
 
if (!library(A2TEA.WebApp, logical.return=T)) 
  quit(status=10)

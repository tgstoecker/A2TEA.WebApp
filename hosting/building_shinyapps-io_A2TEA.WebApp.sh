#!/bin/bash

#make sure that all libraries needed to run the Shiny App are installed in your current (e.g. dev. environment) environment
#alternatively, activate conda environment with all A2TEA.WebApp libraries installed + rsconnect connected to shinyapps.io account
# prior to starting the script

#Step 1 - copying all files shinyapps dir

cp ../inst/webapp/ui.R .
cp ../inst/webapp/server.R .
cp ../R/global.R .
cp ../inst/webapp/www/a2tea_hexsticker.png www/
cp ../inst/webapp/test_data/test.RData .

#Step 2 - modify files to a format suitable for shinyapps.io
##remove all import statements and add 3 library calls
###remove all comments, class objects & faulty logo path from global.R
grep -v '^#' global.R > tmp && rm global.R && mv tmp global.R
sed -i 's/expanded_OG <- //' global.R
sed -i 's/hypothesis <- //' global.R
sed -i 's/add_OG_set <- //' global.R
sed -i 's/webapp\///' global.R

###remove ui.R - wrong path
sed -i 's/A2TEA.WebApp\///' ui.R

###remove server.R - wrong path
sed -i 's/www\/test.RData/test\.RData/' server.R

#Step 3 - start R (environment with all libs of the App installed!) and upload
Rscript -e 'library(shiny)
            options(repos = BiocManager::repositories())
            rsconnect::deployApp(appName="A2TEA-WebApp")'

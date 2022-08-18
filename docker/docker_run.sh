#!/bin/bash

#Building the A2TEA.WebApp docker image
##We make use of some code from https://github.com/rocker-org/rocker-versioned2, namely the rstudio base image.
##This comes preinstalled with some useful libraries such as remotes.

#build the Dockerfile - clean with latest version of base image
docker build --no-cache . -t tgstoecker/a2tea_webapp

#push the image to dockerhub
docker push tgstoecker/a2tea_webapp:latest

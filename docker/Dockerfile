FROM rocker/rstudio:4.1.3

RUN mkdir -p /scripts
COPY installation.sh /scripts/installation.sh
COPY install2.r /scripts/install2.r
COPY install_a2tea_webapp.r scripts/install_a2tea_webapp.r

WORKDIR /scripts
RUN chmod +x installation.sh
RUN chmod +x install2.r
RUN chmod +x install_a2tea_webapp.r

RUN ./installation.sh

RUN chmod a+rwx -R /usr/local/lib/R/site-library/A2TEA.WebApp/webapp/

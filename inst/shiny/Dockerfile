# DOCKERFILE for GUPPI
# github.com/davidsbutcher/guppi
# docker build --tag davidsbutcher/guppi:releaseX --no-cache .

# get shiny plus tidyverse packages image
FROM rocker/shiny-verse:4.0.3

# Install libraries needed for imagemagick
RUN sudo apt update
RUN sudo apt -y install libmagick++-6.q16-dev

# install R packages from Bioconductor

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install('GO.db')"
RUN R -e "BiocManager::install('UniProt.ws')"

# install R packages from github

RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"

RUN R -e "devtools::install_github('davidsbutcher/GUPPI')"
RUN R -e "devtools::install_github('davidsbutcher/viztools')"
RUN R -e "devtools::install_github('hrbrmstr/waffle')"
RUN R -e "devtools::install_github('rstudio/fontawesome')"

# install R packages from CRAN 

RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magick', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('assertthat', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tictoc', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rcpp', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyWidgets', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyFiles', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyjs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('colourpicker', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('purrr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Peptides', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('forcats', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rmarkdown', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('flexdashboard', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('tibble', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('knitr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('UpSetR', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pander', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DT', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggthemes', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('viridis', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RSQLite', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DBI', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('fs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dbplyr', repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('tidyr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('writexl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('sessioninfo', repos='http://cran.rstudio.com/')"

# copy the app to the image
COPY ui.R /var/guppi/
COPY server.R /var/guppi/
COPY about.md /var/guppi/
COPY /www/ /var/guppi/www/

# select port
EXPOSE 3838

# run app
CMD R -e 'shiny::runApp("/var/guppi", port = 3838, host = "0.0.0.0")'
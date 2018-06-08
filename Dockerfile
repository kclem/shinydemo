FROM rocker/shiny:latest
MAINTAINER Kendell Clement "kclement@mgh.harvard.edu"

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev

COPY . /srv/shiny-server/.

RUN Rscript -e "install.packages( \
    c('httr', 'plotly', 'RColorBrewer'),quiet = TRUE)"


COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /app /srv/shiny-server/


EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh
#ENTRYPOINT "Rscript runGitHub('https://github.com/kclem/shinydemo.git',port=3838)"

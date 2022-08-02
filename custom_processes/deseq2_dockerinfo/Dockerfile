FROM r-base:latest

# change the user for install
USER root

# install python3
RUN apt-get update && apt-get install -y build-essential wget unzip default-jre python3 tabix procps && apt-get clean && apt-get purge

# install dependencies
RUN apt-get install libxml2-dev -y
RUN apt-get install libcurl4-openssl-dev -y
RUN apt-get install libssl-dev -y

# helps clean up Docker images
RUN rm -rf /var/lib/apt/lists/*

# get R installation file
ADD install.R /tmp/

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

# install necessary R packages
RUN R -f /tmp/install.R

# copy over necessary scripts
COPY ./bin/* /usr/local/bin

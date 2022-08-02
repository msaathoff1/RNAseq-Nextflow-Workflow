FROM r-base:latest

USER root

RUN apt-get update && apt-get install -y build-essential wget unzip default-jre python3 tabix procps && apt-get clean && apt-get purge

RUN apt-get install libxml2-dev -y
RUN apt-get install libcurl4-openssl-dev -y
RUN apt-get install libssl-dev -y

# Helps clean up Docker images
RUN rm -rf /var/lib/apt/lists/*

ADD install.R /tmp/

# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/

RUN R -f /tmp/install.R

COPY ./bin/* /usr/local/bin

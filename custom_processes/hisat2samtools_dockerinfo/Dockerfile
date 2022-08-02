FROM biocontainers/samtools:v1.9-4-deb_cv1

# set environment variables
ENV hisat2_version 2.1.0

USER root

# Install dependencies
RUN apt-get update && apt-get install -y build-essential wget unzip default-jre python3 python tabix procps && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# download software
WORKDIR /usr/local/
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v${hisat2_version}.tar.gz

# unpack
RUN tar -xvzf v${hisat2_version}.tar.gz

# # install hisat2
WORKDIR /usr/local/hisat2-${hisat2_version}
RUN make
RUN ln -s /usr/local/hisat2-${hisat2_version}/hisat2 /usr/local/bin/hisat2
RUN ln -s /usr/local/hisat2-${hisat2_version}/hisat2-build /usr/local/bin/hisat2-build

USER biodocker

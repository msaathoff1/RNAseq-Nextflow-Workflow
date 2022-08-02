FROM python:3

# change user
USER root

# copy necessary scripts
COPY ./bin/* /usr/local/bin

# install gffutils python library
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get_pip.py
RUN chmod +x /tmp/get_pip.py && python /tmp/get_pip.py
RUN pip install gffutils

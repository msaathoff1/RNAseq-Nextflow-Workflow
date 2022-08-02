FROM python:3

USER root
COPY ./bin/* /usr/local/bin

ADD https://bootstrap.pypa.io/get-pip.py /tmp/get_pip.py
RUN chmod +x /tmp/get_pip.py && python /tmp/get_pip.py
RUN pip install gffutils

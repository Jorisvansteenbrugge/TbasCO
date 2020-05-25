FROM ubuntu:18.04

COPY . /package_source

RUN apt update

RUN apt install -y  git build-essential wget autoconf

RUN wget https://www.freestatistics.org/cran/bin/linux/ubuntu/bionic-cran35/r-base-core_3.6.3-1bionic_amd64.deb -O /tmp/r-base.deb
RUN wget https://www.freestatistics.org/cran/bin/linux/ubuntu/bionic-cran35/r-base-dev_3.6.3-1bionic_all.deb -O /tmp/r-base-dev.deb

RUN dpkg -i /tmp/r-base.deb
RUN apt install -f -y
RUN dpkg -i /tmp/r-base-dev.deb
RUN apt install -f -y
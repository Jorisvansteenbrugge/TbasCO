From ubuntu:18.04

COPY . /package_source

RUN apt update

RUN apt install -y  git build-essential wget autoconf

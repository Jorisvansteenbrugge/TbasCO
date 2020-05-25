FROM r-base:3.6.3

COPY . /package_source
RUN apt update
RUN apt install -y libcurl4-openssl-dev libxml2-dev libgit2-dev libssl-dev

RUN R -e 'install.packages("devtools")'

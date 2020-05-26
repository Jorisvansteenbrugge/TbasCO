FROM r-base:3.6.3

COPY . /package_source
RUN apt update
RUN apt install -y libcurl4-openssl-dev libxml2-dev libgit2-dev libssl-dev libv8-dev

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("V8")'
RUN R -e 'devtools::install_local("/package_source/")'

FROM rocker/tidyverse:4.1.3

WORKDIR /app
COPY . .

RUN apt-get update && \
    apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-pip

RUN pip install cget

RUN Rscript extdata/install_packages.R

# RUN R CMD INSTALL .
ENTRYPOINT ["/bin/bash"]

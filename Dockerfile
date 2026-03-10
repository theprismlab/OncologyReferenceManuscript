FROM rocker/r-ver:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    cmake \
    pkg-config \
    libfontconfig1-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY=renv/library
RUN Rscript -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "renv::restore()"

COPY scripts/ /app/scripts/

ENTRYPOINT ["Rscript"]

CMD ["scripts/1 - DATA_PROCESSING.R"]
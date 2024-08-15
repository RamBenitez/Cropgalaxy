# Use the official R base image
FROM rocker/r-ver:4.0.2

# Set environment variables for R
ENV R_LIBS=/usr/local/lib/R/site-library
ENV PATH="/usr/local/bin:${PATH}"

# Install required packages
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and required R packages
RUN R -e "install.packages('optparse', repos='http://cran.rstudio.com/')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.rstudio.com/'); BiocManager::install(c('fgsea'))"

# Create app directory
WORKDIR /usr/local/src

# Copy fgsea.R and the rest of the Galaxy tool scripts to the container
COPY fgsea.R /usr/local/src/
COPY fgsea_rice.xml /usr/local/src/

ENTRYPOINT ["Rscript", "/usr/local/src/fgsea.R"]
CMD ["--help"]  # Default command if none is provided

VOLUME ["/data/input", "/data/output"]

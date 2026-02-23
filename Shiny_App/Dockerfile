FROM rocker/shiny

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo libcurl4-gnutls-dev libcairo2-dev libxt-dev libssl-dev \
    libssh2-1-dev libglpk40 libproj-dev && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("DT", "data.table", "dplyr", "tidyr", \ 
    "edgeR", "pheatmap", "ggplot2", "cowplot", "corrplot", \ 
    "plotly", "stringr", "BiocManager", "devtools"), dependencies=TRUE)'

RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("edgeR")'
RUN R -e 'devtools::install_github("kevinblighe/PCAtools")'
RUN R -e 'BiocManager::install("AnnotationDbi")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'

RUN	rm -rf /srv/shiny-server/*

WORKDIR	/srv/shiny-server/

COPY ./ ./

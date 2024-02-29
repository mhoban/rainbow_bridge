FROM bioconductor/r-ver:RELEASE_3_17

MAINTAINER "Mykle Hoban" mykle.hoban@hawaii.edu


RUN /rocker_scripts/install_tidyverse.sh

# install packages
RUN Rscript -e 'install.packages(c("insect","optparse","phangorn","furrr","docopt","worrms","rentrez"))' -e 'devtools::install_github(c("tobiasgf/lulu","GuillemSalazar/EcolUtils"))'

RUN Rscript -e 'BiocManager::install(c("phyloseq","dada2","DECIPHER"))'
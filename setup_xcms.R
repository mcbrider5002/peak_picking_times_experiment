repos = 'https://www.stats.bris.ac.uk/R/'

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos=repos)

BiocManager::install("xcms")

install.packages("magrittr", repos=repos)
install.packages("optparse", repos=repos)
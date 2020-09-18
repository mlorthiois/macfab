if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!"bambu" %in% installed.packages()){
    print("Installation of bambu")
    devtools::install_github("GoekeLab/bambu@v0.2.0", quiet=TRUE)
} else {
    print("bambu already installed")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!"BSgenome" %in% installed.packages()){
    "Installation of BSGenome"
    BiocManager::install("BSgenome")
} else {
    print("BSgenome already installed")
}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools", repos="http://cran.irsn.fr/", quiet = TRUE)

if (!requireNamespace("reshape2", quietly = TRUE))
    install.packages("reshape2", repos="http://cran.irsn.fr/", quiet = TRUE)

if (!requireNamespace("tidyverse", quietly = TRUE))
	install.packages("tidyverse", repos="http://cran.irsn.fr/", quiet = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quiet = TRUE)


if (!"bambu" %in% installed.packages()){
    print("Installation of bambu")
    devtools::install_github("GoekeLab/bambu@v0.2.0", quiet=TRUE)
} else {
    print("bambu already installed")
}


if (!"BSgenome" %in% installed.packages()){
    "Installation of BSGenome"
    BiocManager::install("BSgenome")
} else {
    print("BSgenome already installed")
}
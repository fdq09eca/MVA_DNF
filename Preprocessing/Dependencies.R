bioc_packages <- c(
    "PharmacoGx", "annotate",
    "org.Hs.eg.db", "survcomp"
)
carn_packages <- c(
    "BiocManager", "apcluster", "rcdk", "fingerprint",
    "SNFtool", "ROCR", "reshape2",
    "proxy", "pheatmap", "GSA"
)

packages <- c(bioc_packages, carn_packages)

carn_packages <- carn_packages[!(carn_packages %in% installed.packages())]
bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages())]

if (length(carn_packages)) install.packages(carn_packages)
if (length(bioc_packages)) BiocManager::install(bioc_packages)

for (pkg in packages) {
    print(paste("loading package: ", pkg))
    require(pkg, character.only = TRUE)
}

rm(list = ls())
# install packages from CRAN
pkgs.cran <- c(
  "devtools", "remotes", "BiocManager", "ragg", "bit64","lazyeval", "tictoc",
  "argparser", "survival", "zoo", "data.table", "tidyverse",
  "Seurat", "rliger", "harmony", "scCATCH", "SoupX", "NMF"
)
installed_packages <- pkgs.cran %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(pkgs.cran[!installed_packages])
}

# install packages from Bioconductor
pkgs.bioc <- c("ComplexHeatmap","SingleR","clusterProfiler","GSVA",
               "AUCell","DropletUtils", "bambu")
installed_packages <- pkgs.bioc %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(pkgs.bioc[!installed_packages])
}

# install packages from Github
pkgs.github <- c(
  "DoubletFinder", "CellChat", "monocle3", "VISION",
  "scMetabolism", "harmony",
  "SeuratData", "SeuratWrappers", "SeuratDisk",
  "scater", "celldex"
)
pkgs.gitpath <- c(
  "chris-mcginnis-ucsf/DoubletFinder", "sqjin/CellChat",
  "cole-trapnell-lab/monocle3", "YosefLab/VISION",
  "wu-yc/scMetabolism", "immunogenomics/harmony",
  "satijalab/seurat-data", "satijalab/seurat-wrappers",
  "mojaveazure/seurat-disk", "davismcc/scater", "LTLA/celldex"
)
installed_packages <- pkgs.github %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(pkgs.gitpath[!installed_packages])
}


#!/usr/bin/env Rscript
#  Packages installed using bioconductor
pkg_bioc = c("shape","blob","ComplexHeatmap","SingleR","DropletUtils",
             "clusterProfiler","GSVA","AUCell")
#  Packages installed using CRAN
pkg_cran = c("devtools","remotes","ragg","argparser","bit64",
             "survival","zoo","data.table","lazyeval","SoupX",
             "Seurat","rliger","harmony","scCATCH","magick","")
#  Packages installed using github
pkg_gthb = c("DoubletFinder","CellChat","monocle3","VISION","scMetabolism",
             "SCENIC","celldex","SeuratData","SeuratWrappers","SeuratDisk",
             "scater","ggunchull","jjAnno","scRNAtoolVis")
pkg_gthb_install = c("chris-mcginnis-ucsf/DoubletFinder","sqjin/CellChat",
                     "cole-trapnell-lab/monocle3","YosefLab/VISION",
                     "wu-yc/scMetabolism","aertslab/SCENIC","LTLA/celldex",
                     "satijalab/seurat-data","satijalab/seurat-wrappers",
                     "mojaveazure/seurat-disk","davismcc/scater",
                     "sajuukLyu/ggunchull","junjunlab/jjAnno","junjunlab/scRNAtoolVis")

installation = FALSE
tryCatch(
  expr = {
    v = as.numeric(paste0(R.Version()$major, R.Version()$minor)) #36.6

    if(v < 35){
      if (length(setdiff(pkg_cran, rownames(installed.packages()))) > 0) {
        install.packages(setdiff(pkg_cran, rownames(installed.packages())),repos = "https://cloud.r-project.org/")
      }
      if (length(setdiff(pkg_bioc, rownames(installed.packages()))) > 0) {
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(setdiff(pkg_bioc, rownames(installed.packages())), suppressUpdates=TRUE)
      }
      if (length(setdiff(pkg_gthb, rownames(installed.packages()))) > 0) {
        library(devtools)
        devtools::install_github(setdiff(pkg_gthb_install, rownames(installed.packages())))
      }
      installation = TRUE
    }else{
      pkg_cran = c("BiocManager", pkg_cran)
      if (length(setdiff(pkg_cran, rownames(installed.packages()))) > 0) {
        install.packages(setdiff(pkg_cran, rownames(installed.packages())), repos = "https://cloud.r-project.org/")
      }
      if (length(setdiff(pkg_bioc, rownames(installed.packages()))) > 0) {
        lst = setdiff(pkg_bioc, rownames(installed.packages()))
        for(p in lst){
          library(BiocManager)
          BiocManager::install(p, ask=FALSE)
        }
      }
      if (length(setdiff(pkg_gthb, rownames(installed.packages()))) > 0) {
        library(devtools)
        devtools::install_github(setdiff(pkg_gthb_install, rownames(installed.packages())), quiet = T)
      }
    }
    installation = TRUE
  },
  error = function(e){
    installation = FALSE
  }
)
message(installation)

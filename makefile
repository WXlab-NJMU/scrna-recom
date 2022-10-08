SHELL = /bin/bash
OS = $(shell lsb_release -i | cut -d: -f2)
PLATFORM = $(shell uname -s)
RVERSION = $(shell R --version | head -1 | awk '{print $$3}')
# local R library path
RLIB = $(HOME)/miniconda3/envs/scrna/lib/R/library/
all: install test
.PHONY : install clean test


test:
	@if [[ $(strip $(OS)) == "Ubuntu" ]]; then \
		if [[ 1 == 1 ]]; then echo "a"; fi; \
		echo "testing..." ;\
  fi;

conda-r:
$(info "Current R version: $(RVERSION)")
ifeq (,$(shell which R))
	$(info Installing R)
	mamba activate scrna  && mamba install -c conda-forge r-base==4.2 r-essentials r-rcppeigen r-httr pkg-config && conda install -c bioconda bioconductor-rhdf5lib
endif
ifneq (ok, $(shell [[ '$(RVERSION)' > '4.0.0' ]]  && echo ok ))
	$(error "Please update the R version, and it must be greater than 4!")
endif

r-pkgs: conda-r
	#sudo apt-get install libblas-dev liblapack-dev libgeos-dev libcurl4-openssl-dev libhdf5-dev libfontconfig1-dev
	#sudo yum install blas-devel lapack-devel geos-devel libcurl-devel hdf5-devel fontconfig-devel
	$(info Install R Libraries to $(RLIB))
	#mamba install -c conda-forge r-seurat r-seuratdisk r-terra r-stringi r-reshape2
	#mamba install -c bioconda bioconductor-clusterprofiler bioconductor-complexheatmap bioconductor-gsva bioconductor-singler r-monocle3 scvelo r-harmony && mamba install -c bioturing r-rliger r-seuratdata && mamba install -c paul.martin-2 r-doubletfinder
	R --vanilla -e 'install.packages(c("devtools", "remotes", "BiocManager", "ragg", "docopt"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("survival", "zoo", "data.table", "lazyeval", "Seurat", "rliger", "harmony", "scCATCH"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("bit64", "rliger"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("harmony", "scCATCH"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(version = "3.15", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("shape", "blob", "ComplexHeatmap"), lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("SingleR","clusterProfiler","GSVA"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("chris-mcginnis-ucsf/DoubletFinder", "sqjin/CellChat", "wu-yc/scMetabolism", "satijalab/seurat-data"), lib="$(RLIB)")'
	#R --vanilla -e 'remotes::install_github(c("YosefLab/VISION", "aertslab/SCENIC"), lib="$(RLIB)")'
	pip install torch pyscenic
	#pip install cellphonedb


python-pkgs:
ifeq (, $(shell conda))
	$(error Please install conda first!!!)
	exit
endif
ifeq (, $(shell conda env list | grep scrna))
	$(info Conda environment is creating...)
	# mamba create -n scrna python=3.10 && pip install Cython sphinx-book-theme
	mamba env create -f scrna_environment.yaml
endif



SHELL = /bin/bash
OS = $(shell lsb_release -i | cut -d: -f2)
PLATFORM = $(shell uname -s)
RVERSION = $(shell R --version | head -1 | awk '{print $$3}')
# local R library path
RLIB = $(HOME)/anaconda3/envs/scrna/lib/R/library/
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
	conda activate scrna  && conda install -c conda-forge r-base=4.2.1 r-essentials r-biocmanager r-remotes r-rcppeigen r-httr r-docopt r-ragg && conda install -c bioconda bioconductor-rhdf5lib
endif
ifneq (ok, $(shell [[ '$(RVERSION)' > '4.0.0' ]]  && echo ok ))
	$(error "Please update the R version, and it must be greater than 4!")
endif

seurat: conda-r
	#sudo apt-get install libblas-dev liblapack-dev libgeos-dev libcurl4-openssl-dev libhdf5-dev libfontconfig1-dev
	#sudo yum install blas-devel lapack-devel geos-devel libcurl-devel hdf5-devel fontconfig-devel
	$(info Install Seurat to $(RLIB))
	conda activate scrna  && conda install -c conda-forge r-seurat r-seuratdisk r-terra r-stringi r-reshape2
	#conda install -c bioconda bioconductor-clusterprofiler bioconductor-complexheatmap bioconductor-gsva bioconductor-singler r-monocle3 scvelo r-harmony && conda install -c bioturing r-rliger r-seuratdata && conda install -c paul.martin-2 r-doubletfinder && pip install torch pyscenic cellphonedb
	R --vanilla -e 'install.packages(c("Seurat", "rliger", "harmony", "scCATCH"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("ComplexHeatmap","SingleR","clusterProfiler","GSVA"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("chris-mcginnis-ucsf/DoubletFinder", "sqjin/CellChat","YosefLab/VISION", "aertslab/SCENIC", "wu-yc/scMetabolism", "satijalab/seurat-data"), lib="$(RLIB)")'


python:
ifeq (, $(shell conda))
	$(error Please install conda first!!!)
	exit
endif
ifeq (, $(shell conda env list | grep scrna))
	$(info Conda environment is creating...)
	# conda create -n scrna python=3.10 && pip install Cython sphinx-book-theme
	conda env create -f scrna_environment.yaml
endif



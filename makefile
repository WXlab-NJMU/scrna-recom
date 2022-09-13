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
	conda activate scrna && conda install h5py && conda install -c conda-forge r-base=4.2.1 r-essentials r-docopt r-ragg && pip install python-igraph louvain pybind11 hnswlib pyscenic scvelo cellphonedb
endif
ifneq (ok, $(shell [[ '$(RVERSION)' > '4.0.0' ]]  && echo ok ))
	$(error "Please update the R version, and it must be greater than 4!")
endif

seurat: conda-r
	#sudo apt-get install libblas-dev liblapack-dev libgeos-dev libcurl4-openssl-dev libhdf5-dev libfontconfig1-dev
	#sudo yum install blas-devel lapack-devel geos-devel libcurl-devel hdf5-devel fontconfig-devel
	$(info Install Seurat to $(RLIB))
	R --vanilla -e 'install.packages(c("remotes", "BiocManager", "httr", "ploty", "RcppEigen", "Seurat", "rliger", "harmony", "scCATCH"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("HDF5Array", "lme4", "reshape2", "spdep", "stringr", "terra"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("Rhdf5lib","ComplexHeatmap","SingleR","clusterProfiler","GSVA"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("rspatial/terra","mojaveazure/seurat-disk","chris-mcginnis-ucsf/DoubletFinder","sqjin/CellChat","YosefLab/VISION", "aertslab/SCENIC", "wu-yc/scMetabolism", "satijalab/seurat-data","cole-trapnell-lab/monocle3"), lib="$(RLIB)")'


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


clean :
	sudo apt remove r-*
	sudo apt autoremove

rr :
	sudo apt update -y

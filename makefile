SHELL = /bin/bash
OS = $(shell lsb_release -i | cut -d: -f2)
PLATFORM = $(shell uname -s)
RVERSION = $(shell R --version | head -1 | awk '{print $$3}')
# local R library path
RLIB = $(HOME)/miniconda3/envs/scrna-recom/lib/R/library
CRAN = https://mirrors.ustc.edu.cn/CRAN/
all: install
.PHONY : install install-r clean


conda-r:
$(info "Current R version: $(RVERSION)")
ifneq (ok, $(shell [[ '$(RVERSION)' == '4.1.3' ]]  && echo ok ))
	$(info Installing R)
	mamba install -c conda-forge r-base==4.1.3 r-essentials r-rcppeigen r-httr r-gert r-xml r-spdep gcc pkg-config udunits2 gdal lerc && mamba install -c bioconda bioconductor-rhdf5lib
endif

r-pkgs: conda-r
	#sudo apt-get install libblas-dev liblapack-dev libgeos-dev libcurl4-openssl-dev libhdf5-dev libfontconfig1-dev
	#sudo yum install blas-devel lapack-devel geos-devel libcurl-devel hdf5-devel fontconfig-devel
	$(info Install R Libraries to $(RLIB))
	R --vanilla -e 'install.packages(c("devtools", "remotes", "BiocManager", "ragg", "argparser"), repos="$(CRAN)", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("survival", "zoo", "data.table", "lazyeval", "Seurat"), repos="$(CRAN)", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("bit64", "rliger"), repos="$(CRAN)", lib="$(RLIB)")'
	R --vanilla -e 'install.packages(c("harmony", "scCATCH"), repos="$(CRAN)", lib="$(RLIB)")'
	#R --vanilla -e 'BiocManager::install(version = "3.14", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("shape", "blob", "ComplexHeatmap"), lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("SingleR","clusterProfiler","GSVA", "AUCell"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("chris-mcginnis-ucsf/DoubletFinder", "sqjin/CellChat", "cole-trapnell-lab/monocle3"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("YosefLab/VISION", "wu-yc/scMetabolism", "aertslab/SCENIC"), lib="$(RLIB)")'
	R --vanilla -e 'remotes::install_github(c("satijalab/seurat-data", "satijalab/seurat-wrappers", "mojaveazure/seurat-disk", "davismcc/scater", "LTLA/celldex"), lib="$(RLIB)")'
	# install SoupX
	R --vanilla -e 'install.packages(c("SoupX"), repos="$(CRAN)", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("DropletUtils"), lib="$(RLIB)")

python-pkgs:
ifeq (, $(shell conda))
	$(error Please install conda first!!!)
	exit
endif
ifeq (, $(shell conda env list | grep scrna-recom))
	$(info Conda environment is creating...)
	if [ "$(PLATFORM)" == "Linux" ]; then \
		mamba env create -f scrna-recom_linux.yaml ; \
	elif [ "$(PLATFORM)" == "Darwin" ]; then \
		mamba env create -f scrna-recom_macos.yaml ; \
	fi;
endif

install-python: python-pkgs
install-r: r-pkgs


SHELL = /bin/bash
OS = $(shell lsb_release -i | cut -d: -f2)
PLATFORM = $(shell uname -s)
RVERSION = $(shell R --version | head -1 | awk '{print $$3}')
# local R library path
RLIB = $(HOME)/miniconda/envs/scrna/lib/R/library/
all: install test
.PHONY : install clean test


test:
	@if [[ $(strip $(OS)) == "Ubuntu" ]]; then \
		if [[ 1 == 1 ]]; then echo "a"; fi; \
		echo "testing..." ;\
  fi;

R:
# Install R on different OS
	@if [[ "$(shell which R)" == "" ]]; then \
		if [[ "$(strip $(OS))" == "Ubuntu" ]]; then \
			sudo apt update -y; \
			sudo apt install -y --no-install-recommends software-properties-common dirmngr; \
			wget -qo- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && rm marutter_pubkey.asc; \
			sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(shell lsb_release -cs)-cran40/"; \
			sudo apt install -y --no-install-recommends r-base; \
		elif [[ $(PLATFORM) == "Darwin" ]]; then \
			brew install R; \
		fi; \
	fi;
$(info Checking R version...)

conda-R:
ifeq (,$(shell which Rt))
	$(info Installing R)
	conda install -c r r-base=4.2.1 r-essentials r-docopt scvelo cellphonedb
endif
ifneq (ok, $(shell [[ '$(RVERSION)' > '4.0.0' ]]  && echo ok ))
	$(error "Please update the R version, and it must be greater than 4!")
endif

Seurat: R
	sudo apt-get install libblas-dev liblapack-dev libgeos-dev libcurl4-openssl-dev libhdf5-dev libfontconfig1-dev
	$(info Install Seurat to $(RLIB))
	R --vanilla -e 'install.packages(c("remotes", "BiocManager", "httr", "ploty", "RcppEigen", "Seurat", "rliger", "harmony", "scCATCH"), repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'
	R --vanilla -e 'BiocManager::install(c("ComplexHeatmap","SingleR","clusterProfiler","monocle"))'
	R --vanilla -e 'remotes::install_github(c("chris-mcginnis-ucsf/DoubletFinder","sqjin/CellChat","YosefLab/VISION"))'


Python:
ifeq (, $(shell conda))
	$(error Please install conda first!!!)
	exit
endif
ifeq (, $(shell conda env list | grep scrna))
	$(info Conda environment is creating...)
	conda env create -f scrna_environment.yaml
endif


clean :
	sudo apt remove r-*
	sudo apt autoremove

rr :
	sudo apt update -y

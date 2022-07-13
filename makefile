SHELL = /bin/bash
OS = $(shell lsb_release -i | cut -d: -f2)
PLATFORM = $(shell uname -s)
RVERSION = $(shell R --version | head -1 | awk '{print $$3}')
# local R library path
RLIB = $(HOME)/Pipelines/scrna/R-Library/

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
ifneq (ok, $(shell [[ '$(RVERSION)' > '4.0.0' ]]  && echo ok ))
	$(error "Please update the R version, and it must be greater than 4!")
endif

Seurat: R
	$(info Install Seurat to $(RLIB))
	R --vanilla -e 'install.packages("Seurat", repos="https://mirrors.ustc.edu.cn/CRAN/", lib="$(RLIB)")'

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

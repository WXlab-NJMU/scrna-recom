SHELL = /bin/bash
platform := $(shell uname -s)
r_version := $(shell R --version | head -1 | awk '{print $$3}')
install: R
.phony: install clean

R:
ifeq (, $(shell which r))
  $(warning "no r in $(path), r will be installing...")
  ifeq ($(platfrom), linux)
    os := $(shell lsb_release -i | cut -d: -f2)
    os := $(strip $(platform))
    ifeq ($(os), ubuntu)
      $(info "installing r on ubuntu")
      sudo apt update -y
      sudo apt install -y --no-install-recommends software-properties-common dirmngr
      wget -qo- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
      sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
      sudo apt install -y --no-install-recommends r-base
    endif
  endif
  ifeq ($(os), darwin)
    $(info "installing r on mac")
    brew install r
  endif
else
  $(info "Checking R version...")
  ifeq (, $(shell [[ '$(r_version)' > '4.0.0' ]]  && echo ok ))
	  $(info Current R version: $(r_version))
	  $(error "Please update the R version, and it must be greater then 4!")
  endif
endif

Seurat: R
	 R --vanilla -e 'install.packages("Seurat", repos="https://mirrors.ustc.edu.cn/CRAN/")'

clean:
	echo "clean"

test:
  $(info testing)

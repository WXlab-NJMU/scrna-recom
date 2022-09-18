Install
========================================

This project is easily installed through the following command:

.. code-block:: shell

   # for pyscrnarecom in shell
   python setup.py install

   # for scrnaRecom in R console
   library(devtools)
   devtools:install_github("Sue9104/scrna", subdir="scrna-recom")




Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   For CentOS users, we provide a quick installation through makefile.
   We **strongly** recommend these tools should be mannually installed
   for different system environments.

   `make install`


python requirements
'''''''''''''''''''''''''

.. code:: shell

   make python-pkgs

r requirements
'''''''''''''''''''''''''

.. code:: shell

   make r-pkgs


other requirements
'''''''''''''''''''''''''

:bcl2fq: https://emea.support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
:cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

.. important:: Make sure you've added these executed commads before running.


To do
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The docker image of this project is on the way.

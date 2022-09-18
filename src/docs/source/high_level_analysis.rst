High Level Analysis
-----------------------------------

Trajectory Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

scVelo
'''''''''''''''''''''''''

.. parsed-literal::

trajectory-scvelo <indir> <outdir> <h5ad> <gtf> [options]

positional arguments:
  --indir                cellranger directory, including outs and outs/analysis
  --outdir               output directory
  --h5ad                 seurat object h5ad file
  --gtf                  gene gtf file, which is located in ref_folder/genes/genes.gtf

options:
  -h, --help           show this help message and exit
  --genes   gene list, seprate by space
  --mode    stochastic or dynamical
  --start   barcode for start cell


Egulatory Network Inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pySCENIC
'''''''''''''''''''''''''

.. parsed-literal::

network-pysenic <infile> <outdir> <version> <tfs> <cisdb> <motif> [options]

positional arguments:
  --infile               input file
  --outdir               output file
  --version              genome version name
  --tfs                  genome transcript factor file
  --cisdb                cisTarget databse
  --motif                genome motif file

options:
  -h, --help   show this help message and exit
  --minGenes   min gene counts for one cell
  --maxGenes   max gene counts for one cell
  --minCells   min cell counts for one gene
  --maxMT      max percentage of mitochondrial

Metabolic flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

scFEA
'''''''''''''''''''''''''

.. parsed-literal::

metabolism-scfea <infile> <outdir> [options]

positional arguments:
  --infile                Input seurat count matrix csv file
  --outdir                output directory

options:
  -h, --help            show this help message and exit
  --refdir              scFEA model directory
  --species             human or mouse, default is human
  --moduleGene
                        table contains genes for each module
  --stoichiometry       table contains relations between compounds and modules
  --cName               The name of compounds, 2 rows: compounds name and id
  --imputation
                        Whether perform imputation for SC data, default is True

Basic Data Analysis
-----------------------------------

Basic data analysis from count matrix to cluster, including following steps:

- normalization
- feature selection (optional)
- scaling (optional)
- dimensional reduction
- clustering

.. note:: Tools may use different strategies to reduction and clustering, thus steps may be different.

Wrapped tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:R: Seurat https://github.com/satijalab/seurat
:Python: Scanpy https://github.com/scverse/scanpy


Usages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


R: using Seurat in scrnaRecom

.. parsed-literal::
   basic-seruat.R <indir> <outdir> [options]

.. important:: `basic-seruat.R` is in the exec dir of scrnaRecom installed path

Python: using scanpy in pyscrnarecom

.. parsed-literal::
   basic-scanpy <indir> <outdir> [options]

positional arguments:
  --indir                 cellranger filtered_feature_bc_matrix directory
  --outdir                output directory

options:
  -h, --help            show this help message and exit
  --minGenes            min gene counts for one cell
  --maxGenes            max gene counts for one cell
  --minCells            min cell counts for one gene
  --maxMT               max percentage of mitochondrial
  --markerTesting
                        testing method used to find marker genes
  --focusedGenes
                        gene list, separated by space
  --pcNums              run clustering and umap with the fixed pc counts
  --neighbors           neighbors used to find cluster

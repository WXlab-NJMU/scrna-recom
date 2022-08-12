Preprocessing
========================================

Raw Data Processing
-----------------------------------
CellRanger


Quality Control
-----------------------------------
:QC: Seurat https://github.com/satijalab/seurat
:Remove Doublet: DoubletFinder https://github.com/chris-mcginnis-ucsf/DoubletFinder

.. warning::
   Input must be a **fully-processed Seurat object** (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).

   Do not apply DoubletFinder to **aggregated scRNA-seq data** representing multiple distinct samples.

.. code-block:: r
   :linenos:
   :dedent:
   :lineno-start: 0
   :caption: R GUI installation
      install.packages('Seurat')
      remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

Data Integration
-----------------------------------
:Data Integration: Harmony https://github.com/immunogenomics/harmony

.. code-block:: r
   :linenos:
   :dedent:
   :lineno-start: 0
   :caption: R GUI installation
      BiocManager::install("harmony")
      install.packages(c("harmony","rlinger"))




Data Correction
-----------------------------------

Normalization
-----------------------------------
:Normalization: Seurat https://github.com/satijalab/seurat

Feature Selection
-----------------------------------
:Feature Selection: Seurat https://github.com/satijalab/seurat

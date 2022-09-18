Processing
========================================

Raw Data Processing
-----------------------------------
:cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-CellRanger

Convert raw data to count matrix for further analysis
   format: BCL -> Fastq -> Count

   steps: mkfastq -> count -> aggr

.. note:: For cell multiplex barcode, multi is used instead of count.

Requires
'''''''''''''''''''''''''

* bcl2fastq

* cellranger

.. important:: Please make sure these tools are installed and exported to PATH before running!

Usages
'''''''''''''''''''''''''

Here we integrate **multiple cellranger steps** into one command: raw-data-process


.. parsed-literal::

   basic-rawdata <bcl_dir> <outdir> <project> <genome> [options]

-h, --help         show this help message and exit
--dformat          data format: BCL or Fastq
--index            barcode csv, header is Lane,Sample,Index
--fq_path          fastq csv, header is Lane,Sample,FastqDir
--multi            whether to run multi, please set True if using cell multiplex
--aggr             whether to run aggr(normalize,reduction,cluster), please set True if using Loupe Viewer
--aggr_csv         aggr csv, header is sample_id,molecule_h5


Quality Control
-----------------------------------
:QC: Seurat https://github.com/satijalab/seurat
:Remove Doublet: DoubletFinder https://github.com/chris-mcginnis-ucsf/DoubletFinder

.. warning::
   Input must be a **fully-processed Seurat object** (i.e., after NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE have all been run).

   Do not apply DoubletFinder to **aggregated scRNA-seq data** representing multiple distinct samples.


Data Integration
-----------------------------------
:Data Integration: Harmony https://github.com/immunogenomics/harmony
:Data Integration: Liger https://github.com/welch-lab/liger


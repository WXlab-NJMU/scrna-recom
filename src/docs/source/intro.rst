What's scRNA?
-----------------------------------


   A collected resource for scRNA-seq data analysis with biomedical applications

It is challenging for biomedical researchers without bioinformatics background
to understand every detail in scRNA-seq data analysis and conduct data analysis
for their own samples. For instance, scRNA-seq data analysis requires
installation of specific software tools and running through the scripts written
with programming languages such as R and Python.

Along with the recommended workflow, we also provide example computational
scripts together with the software environment setting, which may facilitate
researchers to conduct the data analysis locally.


Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: https://github.com/Sue9104/scrna/blob/main/img/workflow.png?raw=true
   :align: center
   :alt: scRNA recommned workflow

Framework
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tools are integrated in two packages acoording to the program language.

.. list-table:: Framework
   :widths: 5 10 30
   :align: left
   :header-rows: 1
   :stub-columns: 0

   * - Language
     - PackageName
     - WrappedTools
   * - R
     - scrnaRecom
     - Seurat、DoubletFinder、Liger、Harmony、scCATCH、singleR、Monocle3、CellChat、scMetabolism
   * - Python
     - pyscrnarecom
     - cellranger、scanpy、pySCENIC、scVelo、scFEA


Current wrapped tools
'''''''''''''''''''''''''

.. list-table:: Wrapped Tools
   :widths: 15 10 30
   :align: left
   :header-rows: 1
   :stub-columns: 0

   * - Steps
     - Tools
     - Tutorials
   * - Raw Data Processing
     - Cell Ranger
     - https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
   * - Quality Control
     - DoubletFinder
     - https://github.com/chris-mcginnis-ucsf/DoubletFinder
   * -
     - Seurat
     - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
   * - Normalization
     - Seurat
     - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
   * - Integration
     - CCA or RPCA in Seurat
     - https://satijalab.org/seurat/articles/integration_rpca.html
   * -
     - Liger
     - https://github.com/welch-lab/liger
   * -
     - Harmony
     - https://github.com/immunogenomics/harmony
   * - Dimensional Reduction
     - PCA and UMAP in Seurat
     - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#perform-linear-dimensional-reduction-1
   * - Clustering
     - Seurat
     - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells-1
   * - Cell type annotation
     - SingleR
     - https://github.com/dviraran/SingleR
   * -
     - scCATCH
     - https://github.com/ZJUFanLab/scCATCH
   * - Regulon analysis
     - SCENIC
     - https://pyscenic.readthedocs.io/en/latest/
   * - Trajectory inference
     - Monocle3
     - https://cole-trapnell-lab.github.io/monocle3/
   * -
     - scVelo
     - https://github.com/theislab/run_scvelo
   * - Cell communication
     - CellChat
     - https://github.com/sqjin/CellChat
   * - Metabolic flux estimation
     - scMetabolism
     - https://github.com/wu-yc/scMetabolism
   * -
     - scFEA
     - https://github.com/changwn/scFEA


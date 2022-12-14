# scrna-recom

> A collected resource for scRNA-seq data analysis with biomedical applications 

It is challenging for biomedical researchers without bioinformatics background 
to understand every detail in scRNA-seq data analysis and conduct data analysis 
for their own samples. For instance, scRNA-seq data analysis requires 
installation of specific software tools and running through the scripts written 
with programming languages such as R and Python.


Along with the recommended workflow, we also provide example computational 
scripts together with the software environment setting, which may facilitate 
researchers to conduct the data analysis locally. 

**Instructions with practical examples** can be found at:
-  https://github.com/WXlab-NJMU/scrna-recom/wiki

**Complete list of tools** in the paper can be found at:
- https://github.com/WXlab-NJMU/scrna-recom/blob/main/misc/tools.md

## Workflow

<img src="https://github.com/WXlab-NJMU/scrna-recom/blob/main/img/workflow.png" width="80%" height="100%">

## Framework

- R packages are wrapped in ***scrnaRecom***:
  - *qc*: DoubletFinder, SoupX, Seurat
  - *integration*: Liger and Harmony
  - *normalization, reduction and cluster*: Seurat
  - *cell annotation*: singR and scCATCH
  - *trajectory prediction*: Monocle3
  - *cell communication*: CellChat
  - *metabolic flux*: scMetabolism

- Python packages and executations are wrapped in ***pyscrnarecom***:
  - *rawdata*: CellRanger
  - *qc, normalization, reduction and cluster*: scanpy
  - *regulon analysis*: pySCENIC
  - *trajectory prediction*: scVelo
  - *metabolic flux*: scFEA

## Current wrapped tools
<table>
    <thead>
        <tr>
            <th></th>
            <th>Package</th>
            <th>Tutorial</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan=1>Raw Data Processing</td>
            <td>Cell Ranger</td>
            <td>https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger</td>
        </tr>
        <tr>
            <td rowspan=3>Quality Control</td>
            <td>DoubletFinder</td>
            <td>https://github.com/chris-mcginnis-ucsf/DoubletFinder</td>
        </tr>
        <tr>
            <td>Seurat</td>
            <td>https://satijalab.org/seurat/articles/pbmc3k_tutorial.html</td>
        </tr>
        <tr>
            <td>SoupX</td>
            <td>https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html</td>
        </tr>
        <tr>
            <td rowspan=1>Normalization</td>
            <td>Seurat</td>
            <td>https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1</td>
        </tr>
        <tr>
            <td rowspan=3>Integration</td>
            <td>Seurat:CCA,RPCA</td>
            <td>https://satijalab.org/seurat/articles/integration_rpca.html</td>
        </tr>
        <tr>
            <td>Liger</td>
            <td>https://github.com/welch-lab/liger</td>
        </tr>
        <tr>
            <td>Harmony</td>
            <td>https://github.com/immunogenomics/harmony</td>
        </tr>
        <tr>
            <td rowspan=1>Dimensional Reduction</td>
            <td>Seurat:PCA,UMAP</td>
            <td>https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#perform-linear-dimensional-reduction-1</td>
        </tr>
        <tr>
            <td rowspan=1>Clustering</td>
            <td>Seurat</td>
            <td>https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells-1</td>
        </tr>
        <tr>
            <td rowspan=2>Cell type annotation</td>
            <td>SingleR</td>
            <td>https://github.com/dviraran/SingleR</td>
        </tr>
        <tr>
            <td>scCATCH</td>
            <td>https://github.com/ZJUFanLab/scCATCH</td>
        </tr>
        <tr>
            <td rowspan=1>Regulon analysis</td>
            <td>SCENIC</td>
            <td>https://pyscenic.readthedocs.io/en/latest/</td>
        </tr>
        <tr>
            <td rowspan=2>Trajectory inference</td>
            <td>Monocle3</td>
            <td>https://cole-trapnell-lab.github.io/monocle3/</td>
        </tr>
        <tr>
            <td>scVelo</td>
            <td>https://github.com/theislab/scvelo</td>
        </tr>
        <tr>
            <td rowspan=1>Cell communication</td>
            <td>CellChat</td>
            <td>https://github.com/sqjin/CellChat</td>
        </tr>
        <tr>
            <td rowspan=2>Metabolic analysis</td>
            <td>scMetabolism</td>
            <td>https://github.com/wu-yc/scMetabolism</td>
        </tr>
        <tr>
            <td>scFEA</td>
            <td>https://github.com/changwn/scFEA</td>
        </tr>
    </tbody>
</table>




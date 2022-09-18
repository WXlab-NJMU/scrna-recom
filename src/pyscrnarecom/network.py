import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
import glob

import umap
from MulticoreTSNE import MulticoreTSNE as TSNE

import json
import zlib
import base64

from . import ProcessingUsingScanpy

class NetworkInterferenceByPyscenic(object):
    """Single-Cell rEgulatory Network Inference and Clustering Analysis

    pyScenic is used to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data.
    online tutorial: <https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb>

    Attributes:
        indir (str): cellranger filtered_feature_bc_matrix directory
        outdir (str): output directory
        version (str): genome version name
        tf (str): genome transcript factor file
        cisdb (str): cirTarget database path
        motif (str): genome motif file
    """
    workers = 10

    def __init__(self, indir = None, outdir = None,
                 minGenes = 200, maxGenes = 2500, minCells = 3, maxMT = 5,
                 tf = None, cisdb = None, motif = None, version = None):
        self.indir = indir
        self.outdir = os.path.abspath(outdir)
        os.makedirs(self.outdir, exist_ok=True)
        # scanpy qc
        self.minGenes = minGenes
        self.maxGenes = maxGenes
        self.minCells = minCells
        self.maxMT = maxMT
        # database
        self.version = version
        ## transcript factor
        self.tf = tf
        ## cisTarget
        self.cisdb = cisdb
        ## motif
        self.motif = motif

    def run(self):
        # quality control using scanpy
        obj = ProcessingUsingScanpy(
            indir = self.indir, outdir = self.outdir,
            minGenes = self.minGenes, maxGenes = self.maxGenes,
            minCells = self.minCells, maxMT = self.maxMT,
        )
        adata = obj.quality_control()
        # convert to pysenic input format
        row_attrs = { "Gene": np.array(adata.var.index) }
        col_attrs = {
            "CellID":  np.array(adata.obs.index) ,
            "nGene": np.array(np.sum(adata.X.transpose()>0 , axis=0)).flatten(),
            "nUMI": np.array(np.sum(adata.X.transpose() , axis=0)).flatten(),}
        infile = f"{self.outdir}/pysenic.input.loom"
        lp.create(infile, adata.X.transpose(), row_attrs, col_attrs)
        # STEP 1: Gene regulatory network inference, and generation of co-expression modules
        adj_file =  f"{self.outdir}/pysenic.adjectory_matrix.csv"
        cmd = f"pyscenic grn {infile} {self.tf} -o {adj_file} --num_workers {workers}"
        subprocess.run(cmd, shell=True)
        # STEP 2-3: Regulon prediction aka cisTarget from CLI
        reg_file = f"{self.outdir}/pysenic.regulon_prediction.csv"
        db_names = ' '.join( glob.glob(f"{self.cisdb}/*feather") )
        cmd = f"pyscenic ctx {adj_file} {db_names} \
                --annotations_fname {self.motif} \
                --expression_mtx_fname {infile} \
                --output {reg_file} \
                --mask_dropouts \
                --num_workers {worker} "
        subprocess.run(cmd, shell=True)
        # STEP 4: Cellular enrichment (aka AUCell) from CLI
        fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
        sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
        for i,x in enumerate(percentiles):
            fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
            ax.text(x=x, y=ax.get_ylim()[1],
                    s=f'{int(x)} ({percentiles.index.values[i]*100}%)',
                    color='red', rotation=30, size='x-small',
                    rotation_mode='anchor' )
        ax.set_xlabel('# of genes')
        ax.set_ylabel('# of cells')
        fig.tight_layout()
        output = f"{self.outdir}/pyscenic.output.loom"
        cmd = f"pyscenic aucell {infile} {reg_file}  --output {output} \
                --num_workers {workers}"
        subprocess.run(cmd, shell=True)
        # Visualization of SCENIC's AUC matrix
        lf = lp.connect( output, mode='r+', validate=False )
        auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
        lf.close()
        runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
        dr_umap = runUmap( auc_mtx )
        umap_csv = f"{self.outdir}/pyscenic.umap.csv"
        pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( f"{self.outdir}/pyscenic_umap.csv")
        tsne = TSNE( n_jobs=20 )
        dr_tsne = tsne.fit_transform( auc_mtx )
        pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( f"{self.outdir}/pyscenic_tsne.csv")
        # Integration
        final_integrated_file = f'{self.outdir}/pyscenic.integrated-output.loom'
        regulons = lf.ra.Regulons
        auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
        regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
        rt = meta['regulonThresholds']
        for i,x in enumerate(rt):
            tmp = x.get('regulon').replace("(","_(")
            x.update( {'regulon': tmp} )
        ## Concatenate embeddings (tSNE, UMAP, etc.)
        tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])
        Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
        Embeddings_X = pd.concat( [
                pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
                pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
                dr_tsne['X'] ,
                dr_umap['X']
            ], sort=False, axis=1, join='outer' )
        Embeddings_X.columns = ['1','2','3','4']

        Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
        Embeddings_Y = pd.concat( [
                pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
                pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
                dr_tsne['Y'] ,
                dr_umap['Y']
            ], sort=False, axis=1, join='outer' )
        Embeddings_Y.columns = ['1','2','3','4']
        ## SCENIC regulon thresholds:
        metaJson = load_metajson()
        metaJson["regulonThresholds"] = rt
        for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
            clustDict = {}
            clustDict['id'] = i
            clustDict['description'] = f'Unannotated Cluster {i + 1}'
            metaJson['clusterings'][0]['clusters'].append(clustDict)
        clusterings = pd.DataFrame()
        clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)
        # Assemble loom file row and column attributes
        col_attrs = {
            "CellID": np.array(adata.obs.index),
            "nUMI": np.array(adata.obs['n_counts'].values),
            "nGene": np.array(adata.obs['n_genes'].values),
            "Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
            #"Genotype": np.array(adata.obs['Genotype'].values),
            #"Timepoint": np.array(adata.obs['Timepoint'].values),
            #"Sample": np.array(adata.obs['Sample'].values),
            "Percent_mito": np.array(adata.obs['percent_mito'].values),
            "Embedding": dfToNamedMatrix(tsneDF),
            "Embeddings_X": dfToNamedMatrix(Embeddings_X),
            "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
            "RegulonsAUC": dfToNamedMatrix(auc_mtx),
            "Clusterings": dfToNamedMatrix(clusterings),
            "ClusterID": np.array(adata.obs['louvain'].values)
        }
        row_attrs = {
            "Gene": lf.ra.Gene,
            "Regulons": regulons,
        }
        attrs = {
            "title": "sampleTitle",
            "MetaData": json.dumps(metaJson),
            "Genome": self.version,
            "SCopeTreeL1": "",
            "SCopeTreeL2": "",
            "SCopeTreeL3": ""
        }
        # compress the metadata field:
        attrs['MetaData'] = base64.b64encode(
            zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')
        lp.create(
            filename = final_integrated_file,
            layers=lf[:,:],
            row_attrs=row_attrs,
            col_attrs=col_attrs,
            file_attrs=attrs
        )
        lf.close()

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

def run_pyscenic_cli():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run pySCENIC for Egulatory Network Inference')
    parser.add_argument( 'infile', type=str, help='input  file')
    parser.add_argument( 'outdir', type=str, help='output file')
    parser.add_argument( 'version', type=str, help='genome version name')
    parser.add_argument( 'tfs', type=str, help='genome transcript factor file')
    parser.add_argument( 'cisdb', type=str, help='cisTarget databse')
    parser.add_argument( 'motif', type=str, help='genome motif file')
    parser.add_argument( '--minGenes', type=int, default = 200,
                        help='min gene counts for one cell')
    parser.add_argument( '--maxGenes', type=int, default = 2500,
                        help='max gene counts for one cell')
    parser.add_argument( '--minCells', type=int, default = 3,
                        help='min cell counts for one gene')
    parser.add_argument( '--maxMT', type=int, default = 5,
                        help='max percentage of mitochondrial')
    args = parser.parse_args()

def load_metajson():
    ### metadata
    metaJson = {}
    metaJson['embeddings'] = [
        {"id": -1, "name": f"Scanpy t-SNE (highly variable genes)"},
        {"id": 1,  "name": f"Scanpy UMAP  (highly variable genes)"},
        {"id": 2,  "name": "Scanpy PC1/PC2"},
        {"id": 3,  "name": "SCENIC AUC t-SNE"},
        {"id": 4,  "name": "SCENIC AUC UMAP"},
    ]
    metaJson["clusterings"] = [{
                "id": 0,
                "group": "Scanpy",
                "name": "Scanpy louvain default resolution",
                "clusters": [],
            }]
    metaJson["metrics"] = [
            {"name": "nUMI"},
            {"name": "nGene"},
            {"name": "Percent_mito"}
    ]
    metaJson["annotations"] = [
        {"name": "Louvain_clusters_Scanpy", "values": list(set( adata.obs['louvain'].astype(np.str) ))},
        #{"name": "Genotype", "values": list(set(adata.obs['Genotype'].values))},
        #{"name": "Timepoint","values": list(set(adata.obs['Timepoint'].values))},
        #{"name": "Sample", "values": list(set(adata.obs['Sample'].values))}
    ]
    return metaJson

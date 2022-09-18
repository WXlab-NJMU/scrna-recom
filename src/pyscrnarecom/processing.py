import os
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

def scanpy_default_settings():
    # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.verbosity = 2
    sc.logging.print_versions()
    sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)
    sc.settings.autosave = False
    sc.settings.autoshow = False


class ProcessingUsingScanpy(object):
    """Standard Processing using Scanpy.

    Online Tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

    Attributes:
        indir (str): cellranger filtered_feature_bc_matrix directory
        outdir (str): output directory
        minGenes (int): min gene counts for one cell, default is 200
        maxGenes (int): max gene counts for one cell, default is 2500
        minCells (int): min cell counts for one gene, default is 3
        maxMT (int): max percentage of mitochondrial, default is 5
        markerTesting(str): method to find markers, wilcoxon|t-test|logreg
        pcNums (int): the fixed pc numbers
    """

    def __init__(self, indir = None, outdir = None,
                 minGenes = 200, maxGenes = 2500, minCells = 3, maxMT = 5,
                 neighbors = 10, pcNums = None,
                 markerTesting = 'wilcoxon', focusedGenes = None):
        self.indir = indir
        self.outdir = os.path.abspath(outdir)
        os.makedirs(self.outdir, exist_ok=True)
        self.minGenes = minGenes
        self.maxGenes = maxGenes
        self.minCells = minCells
        self.maxMT = maxMT
        self.markerTesting = markerTesting
        self.focusedGenes = focusedGenes
        self.pcNums = pcNums
        self.neighbors = neighbors

    def quality_control(self):
        """quality control"""
        scanpy_default_settings()
        sc.settings.figdir = self.outdir
        # read raw data
        adata = sc.read_10x_mtx(self.indir,var_names='gene_symbols',cache=False)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                   log1p=False, inplace=True)
        sc.pl.violin(adata, save = ".qc_before.pdf",
                     keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                     jitter=0.4, multi_panel=True)
        # filtering
        sc.pp.filter_cells(adata, min_genes = self.minGenes)
        sc.pp.filter_genes(adata, min_cells=self.minCells)
        adata = adata[adata.obs.n_genes_by_counts < self.maxGenes, :]
        adata = adata[adata.obs.pct_counts_mt < self.maxMT, :]
        sc.pl.violin(adata, save = ".qc_after.pdf", show=False,
                     keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                     jitter=0.4, multi_panel=True)
        adata.write(f"{self.outdir}/scanpy.qc_done.h5ad", compression='gzip')
        pd.DataFrame(
            ann.X.toarray(),index=ann.obs.index, columns=ann.var.index
        ).to_csv(f"{self.outdir}/scanpy_count_matrix.qc_done.csv")
        return adata

    def processing_ended_in_pca(self):
        scanpy_default_settings()
        sc.settings.figdir = self.outdir
        prefix = f"{self.outdir}/scanpy_"
        adata = self.quality_control()
        # normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # identify variable genes and filtering
        sc.pp.highly_variable_genes(adata,
                                    min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata,
                                    save=".highly_variable_genes.pdf")
        adata = adata[:, adata.var.highly_variable]
        # scale
        sc.pp.scale(adata, max_value=10)
        # save size
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        # pca and get PCs for louvain and tsne
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca(adata, save=".raw.pdf")
        sc.pl.pca_variance_ratio(adata, log=True, save=".PCs.pdf")
        adata.write(f"{self.outdir}/scanpy.pca_done.h5ad", compression='gzip')

    def run(self):
        pca_file = f"{self.outdir}/scanpy.pca_done.h5ad"
        if not os.path.exists(pca_file):
            self.processing_ended_in_pca()
        if self.pcNums:
            self.processing_in_fixed_pc(pca_file)
            final_file = f"{self.outdir}/scanpy.pca_done.h5ad"
            # plot focused genes
            if self.focusedGenes:
                plot_genes(outdir = self.outdir, h5ad = final_file,
                           genes = self.focusedGenes, keyword = "focused")

    def processing_in_fixed_pc(self, h5ad):
        scanpy_default_settings()
        sc.settings.figdir = self.outdir
        # read
        adata = sc.read(h5ad)
        adata.uns['log1p']["base"] = None
        sc.pp.neighbors(adata, n_neighbors=self.neighbors, n_pcs=self.pcNums)
        # umap
        sc.tl.umap(adata)
        sc.pl.umap(adata, save=f".pc={self.pcNums}.pdf")
        sc.tl.leiden(adata)
        # find marker genes
        ## method: t-test, wilcoxon, logreg
        sc.tl.rank_genes_groups(adata, 'leiden', method=self.markerTesting)
        sc.pl.rank_genes_groups(
            adata, n_genes=25, sharey=False,
            save = ".marker_genes_" + self.markerTesting + ".pdf")
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        df = pd.DataFrame({group + '_' + key[:1]: result[key][group]
                           for group in groups for key in ['names', 'pvals']})
        df.to_csv(f"{self.outdir}/scanpy.find_marker_genes.detail.csv")
        adata.write(f"{self.outdir}/scanpy.final.h5ad", compression='gzip')


def plot_genes(h5ad, outdir, genes, keyword):
    """Plot Focused Genes among clusters"""
    adata = sc.read(h5ad)
    os.makedirs(outdir, exist_ok=True)
    # scanpy settings
    scanpy_default_settings()
    sc.settings.figdir = outdir
    # plot details
    sc.pl.pca(adata, color=genes, save = f".{keyword}.pca.pdf")
    sc.pl.umap(adata, color = ['leiden'] + genes,
               use_raw=False, save = f".{keyword}.umap.pdf")
    sc.pl.violin(adata, genes, groupby='leiden',
                 save = f".{keyword}.dist_in_clusters.violin.pdf")
    sc.pl.dotplot(adata, genes, groupby='leiden',
                 save = f".{keyword}.dist_in_clusters.dotplot.pdf")
    sc.pl.stacked_violin(
        adata, genes, groupby='leiden', rotation = 90,
        save = f".{keyword}.dist_in_clusters.stacked_violin.pdf")

def plot_genes_cli():
    import argparse
    parser = argparse.ArgumentParser( description='Plot Genes after Scanpy Finished')
    group.add_argument('h5ad', type=str, help='h5ad file finished scanpy analysis')
    parser.add_argument( 'outdir', type=str, help='output file')
    parser.add_argument( '--genes', nargs='*', help='gene list, separated by space')
    parser.add_argument( '--keyword', type=str, default='focused',
                        help='output file prefix')
    args = parser.parse_args()
    plot_genes(args.outdir, args.genes, args.keyword, h5ad=args.h5ad)

def run_scanpy_cli():
    import argparse
    parser = argparse.ArgumentParser( description='Run Scanpy Standard Analysis')
    parser.add_argument( 'indir', type=str,
                        help='cellranger filtered_feature_bc_matrix directory')
    parser.add_argument( 'outdir', type=str, help='output directory')
    parser.add_argument( '--minGenes', type=int, default = 200,
                        help='min gene counts for one cell')
    parser.add_argument( '--maxGenes', type=int, default = 2500,
                        help='max gene counts for one cell')
    parser.add_argument( '--minCells', type=int, default = 3,
                        help='min cell counts for one gene')
    parser.add_argument( '--maxMT', type=int, default = 5,
                        help='max percentage of mitochondrial')
    parser.add_argument( '--markerTesting', type=str, default = 'wilcoxon',
                        help='testing method used to find marker genes')
    parser.add_argument( '--focusedGenes', nargs='*', default = None,
                        help='gene list, separated by space')
    parser.add_argument( '--pcNums', type=int, default = None,
                        help='run clustering and umap with the fixed pc counts')
    parser.add_argument( '--neighbors', type=int, default = 10,
                        help='neighbors used to find cluster')
    args = parser.parse_args()
    args = vars(args)
    print(args)

    obj = ProcessingUsingScanpy(**args)
    obj.run()

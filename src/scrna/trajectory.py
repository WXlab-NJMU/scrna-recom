import os
import subprocess
import glob
import logging
import scvelo as scv
class RunSCVeloNow(object):
    """Running scVelo from seurat input

    scvelo online tutorial: https://scvelo.readthedocs.io/

    Attributes:
        h5ad (str): seurat h5ad file
        loom (str): loom file from velocyto pipeline
        outdir (str): output directory
        genes (list): list of focused genes, default is None
        mode (str): stochastic or dynamical, default is stochastic
        start (str): barcode for start cell, default is None
    """
    def __init__(self, h5ad, loom, outdir,
                 genes = None, mode = 'stochastic', start = None):
        self.h5ad = h5ad
        self.loom = loom
        self.outdir = os.path.abspath(outdir)
        os.makedirs(self.outdir, exist_ok=True)
        self.genes = genes
        self.mode = mode
        self.start = start

    def run(self):
        """Run Trajectory Analysis..."""
        scv.logging.print_version()
        scv.settings.verbosity = 3
        scv.settings.presenter_view = False
        scv.set_figure_params('scvelo')
        adata = scv.read(self.h5ad)
        ldata = scv.read(self.loom)
        adata = scv.utils.merge(adata, ladata)
        # set cluster to seurat_clusters
        cluster = 'seurat_clusters'
        if mode == "stochastic":
            self.run_in_stochastic_mode(adata)
        if mode == "dynamical":
            self.run_in_dynamical_mode(adata)

    def run_in_stochastic_mode(self, adata):
        """Run scVelo in stochastic mode"""
        prefix = os.path.abspath(self.outdir) + "/scvelo_stochastic_"
        scv.pl.proportions(adata, groupby = cluster,
                           save = prefix + "proportions.pdf")
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

        scv.tl.velocity(adata, mode="stochastic")
        scv.tl.velocity_graph(adata)

        scv.pl.velocity_embedding_stream(adata, basis='umap',
                                         save = prefix + "embedding-stream.pdf")
        scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120,
                                  save = prefix + "embedding-stream.pdf")
        # velocity for focused genes
        if genes is not None:
            scv.pl.velocity(adata, genes, ncols=2,
                            save = prefix + "genes-velocity.pdf")
            scv.pl.scatter(adata, genes, color=[cluster, 'velocity'],
                           save = prefix + "genes-velocity-scatter.pdf")
        # identify import genes
        scv.tl.rank_velocity_genes(adata, groupby = cluster, min_corr=.3)
        df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
        df.to_csv(prefix + "rank-velocity-genes.csv")
        # speed and coherence
        scv.tl.velocity_confidence(adata)
        keys = 'velocity_length', 'velocity_confidence'
        scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],
                       save = prefix + "speed-coherence.pdf")
        df = adata.obs.groupby(cluster)[keys].mean().T
        df.style.background_gradient(cmap='coolwarm', axis=1).to_excel(
            prefix + "speed-coherence.xls"
        )
        # pseudotime
        ## start is cellbarcode or adata index
        if start is not None:
            x, y = scv.utils.get_cell_transitions(adata, basis='umap',
                                                  starting_cell=start)
            ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05,
                                       show=False)
            ax = scv.pl.scatter(adata, x=x, y=y, s=120, ax=ax,
                                c='ascending', cmap='gnuplot',
                                save = prefix + "cell-pseduotime.pdf")
        scv.tl.velocity_pseudotime(adata)
        scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',
                       save = prefix + "pseudotime.pdf")
        # PAGA velocity
        adata.uns['neighbors']['distances'] = adata.obsp['distances']
        adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
        scv.tl.paga(adata, groups=cluster, save = prefix + "paga-velocity.pdf")
        df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
        df.style.background_gradient(cmap='Blues').format('{:.2g}')
        df.to_excel(prefix + "paga-transitions-confidence.xls")


    def run_in_dynamical_mode(self, adata):
        """Run scVelo in stochastic mode"""
        prefix = os.path.abspath(self.outdir) + "/scvelo_dynamical_"
        scv.tl.velocity(adata, mode='dynamical')
        scv.tl.velocity_graph(adata)
        adata.write(prefix + "dynamical-data.h5ad", compression='gzip')
        # reload
        # adata = scv.read(prefix + "dynamical-data.h5ad")
        scv.pl.velocity_embedding_stream(
            adata, basis='umap',
            save = prefix + "embedding.dynamical.pdf")
        # kinetic rate
        df = adata.var
        df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
        kwargs = dict(xscale='log', fontsize=16)
        with scv.GridSpec(ncols=3) as pl:
            pl.hist(df['fit_alpha'], xlabel='transcription rate',
                    **kwargs, pdf = prefix + "dynamical-transcription-rate.pdf")
            pl.hist(df['fit_beta'] * df['fit_scaling'],
                    xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs,
                    pdf = prefix + "dynamical-splicing-rate.pdf")
            pl.hist(df['fit_gamma'], xlabel='degradation rate',
                    xticks=[.1, .4, 1], **kwargs,
                    pdf = prefix + "dynamical-degradation-rate.pdf")
        scv.get_df(adata, 'fit*', dropna=True).to_csv(
            prefix + "dynamical-kinetic-rates.csv")
        scv.tl.latent_time(adata)
        scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,
                       save = prefix + "latent-time.pdf")
        # top-likelyhood genes
        top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
        scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False,
                       save = prefix + "dynamical-top-liklyhood-genes.pdf")
        scv.tl.rank_dynamical_genes(adata, groupby=cluster)
        df = scv.get_df(adata, 'rank_dynamical_genes/names')
        df.to_csv(prefix + 'rank_dynamical_genes.csv')


def run_velocyto(indir, gtf):
    """Run Velocyto to Generate Spliced and Unspliced Data

    The velocyto pipeline takes a long time to complete and requires high amounts of resources (include RAM/CPUs).
    Find more on <https://velocyto.org/velocyto.py/>
    Install velocyto first: `pip install velocyto`

    Args:
        indir (str): input cellranger directory, including outs, outs/analysis ..
        gtf (str): gene gtf file, which is located in cellranger_ref/genes/genes.gtf
        outdir (str): output directory

    Returns:
        loom files
    """
    cmd = f"velocyto run10x {indir} {gtf}"
    subprocess.run(cmd, shell=True)

def run_scvelo():
    """Trajectory Analysis using scVelo

    scVelo trajectory is based on the balance of spliced and unsplice data.

    Args:
        indir (str): input cellranger directory, including outs, outs/analysis ..
        gtf (str): gene gtf file, which is located in cellranger_ref/genes/genes.gtf
        h5ad (str): seurat h5ad file
        outdir (str): output directory
        genes (str): list of focused genes
        mode (str): stochastic or dynamical, default is stochastic
        start (str): barcode for start cell
    """
    import argparse
    parser = argparse.ArgumentParser( description='scVelo Trajectory Analysis')
    parser.add_argument( 'indir', type=str, help='cellranger directory, including outs and outs/analysis')
    parser.add_argument( 'outdir', type=str, help='output directory')
    parser.add_argument( 'h5ad', type=str, help='seurat object h5ad file')
    parser.add_argument( 'gtf', type=str, help='gene gtf file, which is located in ref_folder/genes/genes.gtf')
    parser.add_argument( '--genes', nargs='*', help='gene list, seprate by space')
    parser.add_argument( '--mode', default='stochastic', help='stochastic or dynamical')
    parser.add_argument( '--start', help='barcode for start cell')
    args = parser.parse_args()
    print(args)
    #genes = None, mode = 'stochastic', start = None):
    logging.info("Running Velocyto...")
    run_velocyto(args.indir, args.gtf)
    loom = glob.glob("indir/*loom")[0]
    logging.info("Running Velocyto...")
    RunSCVeloNow(h5ad = args.h5ad, loom = loom, outdir = args.outdir,
                 genes = args.genes, start = args.start, mode = args.mode)

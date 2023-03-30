import scanpy as sc
from scipy import sparse
import seaborn as sns
import matplotlib
import csv


def plot_cluster_markers(infile, markerfile, outdir, cluster="seurat_clusters"):
    """Plot Cluster Markers using scanpy"""
    seu = sc.read_h5ad(infile)
    seu.obs["seurat_clusters"] = seu.obs["seurat_clusters"].astype("str").astype("category")
    marker_dict = {}
    with open(markerfile, newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene = row["gene"]
            if gene in seu.obs.columns: # transcripts
            #if gene in seu.var.index: #genes
                if row["cellType"] in marker_dict.keys():
                    marker_dict[row["cellType"]] += [row["gene"]]
                else:
                    marker_dict[row["cellType"]] = [row["gene"]]
    print(marker_dict)
    # error if markers in both seu.obs and seu.raw.var_names
    for gene in seu.obs:
        if gene in seu.raw.var_names:
            seu.obs.drop(columns = gene, inplace=True)


    sc.settings.figdir = outdir
    sc.tl.dendrogram(seu, 'seurat_clusters')
    sc.pl.dendrogram(seu, 'seurat_clusters', save="_clusters.pdf")
    sc.pl.umap(
        seu, color='seurat_clusters', add_outline=True, frameon=False,
        legend_loc='on data',legend_fontsize=8, legend_fontoutline=1,
        title='clustering', palette='Spectral', save="_clusters.pdf"
    )
    sc.pl.umap(
        seu, color='orig.ident', add_outline=True, frameon=False,
        legend_loc='right margin',legend_fontsize=8, legend_fontoutline=1,
        title='sample clusters', palette='Spectral', save="_clusters-samples.pdf"
    )
    sc.pl.matrixplot(seu, marker_dict, 'seurat_clusters', cmap='Blues',
                     dendrogram=True, save="_clusters-markers.pdf")
    sc.pl.dotplot(seu, marker_dict, 'seurat_clusters', dendrogram=True,
                  cmap='Blues', save="_clusters-markergene.pdf")
    sc.pl.stacked_violin(seu, marker_dict, groupby='seurat_clusters',
                         swap_axes=False, cmap='Blues', dendrogram=True,
                         save="_clusters-markergene.pdf")
    sc.tl.rank_genes_groups(seu, groupby='seurat_clusters', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap="Blues",
                                    save="_clusters-expression.pdf")
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap='bwr',
                                    values_to_plot='logfoldchanges',
                                    min_logfoldchange=3, vmax=7, vmin=-7,
                                    save="_clusters-logfoldfchange.pdf")

    sc.tl.dendrogram(seu, cluster)
    sc.pl.dendrogram(seu, cluster, save="_celltype.pdf")
    sc.pl.umap(
        seu, color=cluster, add_outline=True, frameon=False,
        legend_loc='on data',legend_fontsize=8, legend_fontoutline=1,
        title='clustering', palette='Spectral', save="_celltype.pdf"
    )
    sc.pl.matrixplot(seu, marker_dict, cluster, cmap='Blues',
                     dendrogram=True, save="_celltype-markers.pdf")
    sc.tl.rank_genes_groups(seu, groupby=cluster, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap="Blues",
                                    save="_celltype-expression.pdf")
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap='bwr',
                                    values_to_plot='logfoldchanges',
                                    min_logfoldchange=3, vmax=7, vmin=-7,
                                    save="_celltype-logfoldfchange.pdf")


def plot_clusters(infile, outdir):
    """Plot Cluster Markers using scanpy"""
    seu = sc.read_h5ad(infile)
    seu.obs["seurat_clusters"] = seu.obs["seurat_clusters"].astype("str").astype("category")
    for gene in seu.obs:
        if gene in seu.raw.var_names:
            seu.obs.drop(columns = gene, inplace=True)

    sc.settings.figdir = outdir
    sc.pl.umap(
        seu, color='seurat_clusters', add_outline=True, frameon=False,
        legend_loc='on data',legend_fontsize=8, legend_fontoutline=1,
        title='clustering', palette='Spectral', save="_clusters.pdf"
    )
    sc.pl.umap(
        seu, color='orig.ident', add_outline=True, frameon=False,
        legend_loc='on data',legend_fontsize=8, legend_fontoutline=1,
        title='sample clusters', palette='Spectral', save="_samples.pdf"
    )
    sc.tl.dendrogram(seu, 'seurat_clusters')
    sc.pl.dendrogram(seu, 'seurat_clusters', save="_clusters.pdf")

    sc.tl.rank_genes_groups(seu, groupby='seurat_clusters', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap="Blues",
                                    save="_expression.pdf")
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap="Blues",
                                    values_to_plot="log10_pvals_adj",
                                    save="_padj.pdf")
    sc.pl.rank_genes_groups_dotplot(seu, n_genes=5, cmap='bwr',
                                    values_to_plot='logfoldchanges',
                                    min_logfoldchange=3, vmax=7, vmin=-7,
                                    save="_logfoldfchange.pdf")
    sc.pl.rank_genes_groups_stacked_violin(seu, n_genes=5, cmap='Blues',
                                           save="_expression.pdf")
    sc.pl.rank_genes_groups_heatmap(seu, n_genes=5, use_raw=False,
                                    swap_axes=True, vmin=-3, vmax=3, cmap='bwr',
                                    show_gene_labels=True, figsize=(16,9),
                                    show=False, save="_logfoldchange.pdf")


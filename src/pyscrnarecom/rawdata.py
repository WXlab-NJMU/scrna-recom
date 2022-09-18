import pandas as pd
import subprocess
from shutil import which
import logging

class RawDataProcessing:

    """Convert Raw Data to Count using CellRanger

    Online Manual: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
    Requirements:
    - bcl2fastq: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
    - cellranger: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

    Attributes:
        indir (str): Input directory
        outdir (str): Output directory
        project (str): Project name
        dformat (str): BCL or Fastq
        index (str): Need when data format is BCL, header is Lane,Sample,Index
        fq_path (str): Need fastq path csv, header is Sample,FastqDir
        multi (bool): execute multi instead of count for cell multiplex barcode, default is False
        aggr (bool): Perform normalize, dimensional reduction and clustering for Loupe Viewer

    Raises:
        ModuleNotFoundError Need install bcl2fastq and cellranger first!!!
        OSError Need different files to run!!!
        NotImplementedError CellRanger multi is not implemented!!!
    """
    def __init__(self, indir=None, outdir=None, project=None, genome=None,
                 dformat = "BCL", index = None, fq_path = None,
                 multi = False, aggr = False, aggr_csv = None):
        self.indir = indir
        self.outdir = outdir
        self.project = project
        self.genome = genome
        self.dformat = dformat
        self.index = index
        self.fq_path = fq_path
        self.multi = multi
        self.aggr = aggr
        self.aggr_csv = aggr_csv

    def run(self):
        """Run CellRanger Analysis"""
        # check whether cmd exists
        requires = ['bcl2fastq', 'cellranger']
        for cmd in requires:
            if which(cmd) is None:
                raise ModuleNotFoundError(f'{cmd} NOT FOUND, Please install first!!!')
        # BCL to Fastq
        if self.type == "BCL":
            if not self.index:
                raise OSError('index csv is needed!', self.index)
            else:
                cmd = f"cellranger mkfastq --id={self.project} \
                    --run={self.bcl} --csv={self.index}"
                subprocess.run(cmd, shell=True)
                ## outputs: {self.project}/outs/fastq_path/*/sample/*fq.gz
                dirs = list(set([ os.path.dirname(folder)
                                 for file in glob.glob(f"{self.bcl}/../{self.project}/outs/fastq_path/*/*/*fastq.gz")]))
                fastqs = [{'Sample': folder.split('/')[-1], 'FastqDir': folder} for folder in dirs]
                self.fq_path = f"{self.bcl}/../{self.project}/sample_fastq.csv"
                pd.DataFrame(fastqs).to_csv(self.fq_path, index=False)
        # Fastq to Count
        if not self.multi:
            cmd = f"cellranger multi --id {self.project}_multi  --csv {self.index}"
            print(cmd)
            raise NotImplementedError("cellranger multi is not implemented, please run manually!!!")
        else:
            if not self.fq_path:
                raise OSError('fq_path csv is needed!', self.fq_path)
            data = pd.read_csv(self.fq_path)
            samples = data['Sample'].tolist()
            procs = []
            for sample in samples:
                fastqdir = data.loc[data['Sample'] == sample,'FastqDir'].squeeze()
                cmd = f"cellranger count --id={self.project}_{sample}_count --fastqs={fastqdir} --sample={sample} --transcriptome={self.genome}"
                ## output: id/outs/filtered_feature_bc_matrix
                ## output: id/outs/molecule_info.h5
                proc = subprocess.Popen(cmd, shell=True)
                procs.append(proc)
            ## do parallel, check if it's running
            for proc in procs:
                proc.communicate()
            ## output aggr csv
            self.aggr_csv = f'{self.bcl}/../{self.project}/sample_aggr.csv'
            info = [{'sample_id': sample, 'molecule_h5': f"{self.project}_{sample}_count/outs/molecule_info.h5"} for sample in samples]
            pd.DataFrame(info).to_csv(self.aggr_csv, index=False)
        # Normalize, Dimensional Reduction and Clustering
        if self.aggr:
            cmd = f"cellranger aggr --id={self.project}_aggr --csv={self.aggr_csv}"
            subprocess.run(cmd, shell=True)
            ## outputs: aggr/outs/count/filtered_feature_bc_matrix
            ## outputs: aggr/outs/count/analysis

            ## outputs: aggr/outs/count/analysis/ clustering,pca,tsne,umap,diffexp

def run_cellranger_cli():
    """Running CellRanger Workflow in one command

    Args:
        indir (str): Input directory
        outdir (str): Output directory
        project (str): Project name
        dformat (str): BCL or Fastq, default is BCL
        index (str): Need when data format is BCL, header is Lane,Sample,Index
        fq_path (str): Need fastq path csv, header is Sample,FastqDir
        multi (bool): execute multi instead of count for cell multiplex barcode, default is False
        aggr (bool): Perform normalize, dimensional reduction and clustering for Loupe Viewer
        aggr_csv (str): Need aggr csv, header is sample_id,molecule_h5

    Raises:
        ModuleNotFoundError Need install bcl2fastq and cellranger first!!!
        OSError Need different files to run!!!
        NotImplementedError CellRanger multi is not implemented!!!

    Returns:
        Analysis Outputs
    """
    import argparse
    parser = argparse.ArgumentParser( description='Raw Data Processing Using CellRanger')
    parser.add_argument( 'indir', type=str, help='input directory')
    parser.add_argument( 'outdir', type=str, help='output directory')
    parser.add_argument( 'project', type=str, help='project name')
    parser.add_argument( 'genome', type=str, help='cellranger genome path')
    parser.add_argument( '--dformat', type=str, help='data format: BCL or Fastq',
                        default = "BCL")
    parser.add_argument('--index', type=str, help="barcode csv, header is Lane,Sample,Index")
    parser.add_argument('--fq_path', type=str, help="fastq csv, header is Lane,Sample,FastqDir")
    parser.add_argument('--multi', type=bool, default = False,
                        help="whether to run multi, please set True if using cell multiplex")
    parser.add_argument('--aggr', type=bool, default = False,
                        help="whether to run aggr(normalize,reduction,cluster), please set True if using Loupe Viewer")
    parser.add_argument('--aggr_csv', type=str, help="aggr csv, header is sample_id,molecule_h5")
    args = parser.parse_args()
    args = vars(args)
    print(args)
    obj = RawDataProcessing(**args)
    obj.run()

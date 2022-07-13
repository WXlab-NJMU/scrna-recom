import subprocess
import pandas as pd
from dotenv import load_dotenv
class DataAnalysis10X(indir, outdir, index):

    """Data Analysis for 10X Genomics: BCL to count

    Attributes:
        indir (str): input bcl folder
        outdir (str): output fastq folder
        index (str): sample index csv file, header is 'Lane,Sample,Index'
    Todo:
        * a1
        * a2

    """

    def __init__(self, indir, outdir):
        """TODO: to be defined. """
        self.indir = indir
        self.outdir = outdir
        self.index = index

    @property
    def samples(self):
        samples = pd.read_csv(self.index)["Sample"].to_list()
        return samples
    @property
    def config(self):
        config = dotenv_values("../.env")
        return config
    @property
    def cell_ranger(self):
        return config["CellRanger"]
    @property
    def cell_ranger(self):
        return config["CellRangerGenome"]
    @property
    def agg_csv(self):
        outdir = self.outdir
        return f"{outdir}/celranger_agg.csv"

    def bcl2fq(self):
        """CellRanger Convert BCL to Fastq Files
        :returns: fastq files
        """
        indir = self.indir
        index = self.index
        genome = self.genome
        cell_ranger = self.cell_ranger
        cmd = f"{cell_ranger} mkfastq --csv={index} \
            --run={indir} --output-dir={outdir}"
        print(cmd)

    def fq2count(self):
        genome = self.genome
        cell_ranger = self.cell_ranger
        header = "sample_id,molecule_h5"
        for sample in self.samples:
            cmd = f"{cell_ranger} --transcriptome={genome} \
                --id={sample}_count --fastqs={outdir} --sample={sample}"
            print(cmd)
            text = f"{sample},{outdir}/{sample}_count/{sample}_molecule_info.h5"
            print(text)


    def count2agg(self):
        csv = self.agg_csv
        cmd = f"cellranger aggr --id=out_aggr --csv={csv}"
        print(cmd)




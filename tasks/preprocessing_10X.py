import os
import subprocess
import pandas as pd
from dotenv import dotenv_values
from typing import List

class DataAnalysis10X:
    """Data Analysis for 10X Genomics: BCL to count

    Attributes:
        indir (str): input bcl folder
        outdir (str): output fastq folder
        index (str): sample index csv file, header is 'Lane,Sample,Index'
    Todo:
        * a1
        * a2
    """

    def __init__(self, indir: str, outdir: str, index: str):
        """TODO: to be defined. """
        self.indir = indir
        self.index = index
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        config = dotenv_values()
        self.cell_ranger = config["CellRanger"]
        self.genome = config["CellRangerGenome"]

    @property
    def samples(self) -> List[str]:
        samples = pd.read_csv(self.index)["Sample"].to_list()
        return samples
    @property
    def agg_csv(self) -> str:
        outdir = self.outdir
        return f"{outdir}/celranger_agg.csv"

    def bcl2fq(self) -> None:
        """CellRanger Convert BCL to Fastq Files
        :returns: fastq files
        """
        cmd = f"{self.cell_ranger} mkfastq --csv={self.index} \
            --run={self.indir} --output-dir={self.outdir}"
        print(cmd)

    def fq2count(self) -> None:
        header = "sample_id,molecule_h5"
        for sample in self.samples:
            cmd = f"{self.cell_ranger} --transcriptome={self.genome} \
                --id={sample}_count --fastqs={self.outdir} --sample={sample}"
            print(cmd)
            text = f"{sample},{self.outdir}/{sample}_count/{sample}_molecule_info.h5"
            print(text)

    def count2agg(self) -> None:
        cmd = f"{self.cell_ranger} aggr --id=out_aggr --csv={self.agg_csv}"
        print(cmd)

    def run(self) -> None:
        self.bcl2fq()
        self.fq2count()
        self.count2agg()

csv = "/home/minsu/Pipelines/scrna/tests/cellranger-tiny-bcl-simple-1.2.0.csv"
bcl = "/home/minsu/Pipelines/scrna/tests/cellranger-tiny-bcl-1.2.0"
DataAnalysis10X(indir=bcl, outdir="test-results", index=csv).run()

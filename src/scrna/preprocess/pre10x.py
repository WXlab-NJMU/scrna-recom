import os
import subprocess
import pandas as pd
from typing import List

class DataAnalysis10X:
    """Data Analysis for 10X Genomics: BCL to count

    Attributes:
        indir (str): input bcl folder
        outdir (str): output fastq folder
        index (str): sample index csv file, header is 'Lane,Sample,Index'
        genome (str): cellranger genome path
    """

    def __init__(self, indir: str, outdir: str, index: str, genome:str):
        """TODO: to be defined. """
        self.indir = indir
        self.index = index
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.genome = genome

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
        cmd = f"cellranger mkfastq --csv={self.index} \
            --run={self.indir} --output-dir={self.outdir}"
        print(cmd)

    def fq2count(self) -> None:
        header = "sample_id,molecule_h5"
        for sample in self.samples:
            cmd = f"cellranger --transcriptome={self.genome} \
                --id={sample}_count --fastqs={self.outdir} --sample={sample}"
            print(cmd)
            text = f"{sample},{self.outdir}/{sample}_count/{sample}_molecule_info.h5"
            print(text)

    def count2agg(self) -> None:
        cmd = f"cellranger aggr --id=out_aggr --csv={self.agg_csv}"
        print(cmd)

    def run(self) -> None:
        self.bcl2fq()
        self.fq2count()
        self.count2agg()

def hello_world():
    """Example function with types documented in the docstring.

    function description

    Args:
            param1 (int): The first parameter.
            param2 (str): The second parameter.

    Returns:
            bool: The return value. True for success, False otherwise.

    Yields:
            int: The next number in the range of 0 to  - 1
            works for generator

    Raises:
        AttributeError: The Raises section is a list of all exceptions
        ValueError: If  is equal to .

    Examples:
            >>> print([i for i in example_generator(4)])
            [0,1,2,3]

    """
    def print():
        print("hi")

class MyClass(object):
    """Summary line for a class

    class description

    Attributes:
        attr1 (str): Description
        attr2 (int, optional): Description

    """
    attr1 = 1
    attr2 = 2



csv = "/home/minsu/Pipelines/scrna/tests/cellranger-tiny-bcl-simple-1.2.0.csv"
bcl = "/home/minsu/Pipelines/scrna/tests/cellranger-tiny-bcl-1.2.0"
#DataAnalysis10X(indir=bcl, outdir="test-results", index=csv).run()

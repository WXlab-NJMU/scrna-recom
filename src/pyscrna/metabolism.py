import subprocess
import os
import logging
import pandas as pd
from shutil import which


class MetabolicFluxUsingSCFEA:

    """cell-wise metabolic estimation using scFEA

    Online Manual:
    - https://github.com/changwn/scFEA/blob/master/scFEA_tutorial1.ipynb
    - https://github.com/changwn/scFEA/blob/master/scFEA_tutorial2.ipynb

    Requirements:
    - scFEA.py: add to PATH, https://github.com/changwn/scFEA

    Attributes:
        infile (str): Input seurat count matrix csv file
        outdir (str): Output directory
        refdir (str): directory including scFEA model files
        species (str): human or mouse
        moduleGene (str): table contains genes for each module, default is module_gene_m168.csv for human
        stoichiometry(str): table contains relations between compounds and modules, default is cmMat_171.csv for human
        cName(str): The name of compounds, 2 rows: compounds name and id
        imputation(bool): Whether perform imputation for SC dataset

    Raises:
        ModuleNotFoundError Need install scFEA first!!!
        OSError Need different files to run!!!
        NotImplementedError CellRanger multi is not implemented!!!
    """
    def __init__(self, infile=None, outdir=None, refdir=None, species=None,
                 moduleGene=None, stoichiometry=None, cName=None,
                 imputation=None):
        self.infile = infile
        self.outdir = outdir
        self.refdir = refdir
        self.species = species
        self.moduleGene = "module_gene_m168.csv" if (species == "human") and (not moduleGene) else moduleGene
        self.stoichiometry = "cmMat_171.csv" if (species == "human") and (not stoichiometry) else stoichiometry
        self.cName = cName
        self.imputation = imputation

    def run(self):
        """Run CellRanger Analysis"""
        # check whether cmd and models exist
        requires = ['scFEA.py']
        for cmd in requires:
            if which(cmd) is None:
                raise ModuleNotFoundError(f'{cmd} NOT FOUND, Please install first!!!')
        scfea = which("scFEA.py")
        if not self.refdir:
            logging.info("scFEA.py path: ", scfea)
            self.refdir = os.path.join(os.path.dirname(scfea), "../data")
            logging.info("scFEA model path: ", self.refdir)
            model = os.path.join(self.refdir,self.moduleGene)
            if not os.path.exists(model):
                raise OSError("scFEA model directory need to exists", model)
        indir = os.path.dirname(os.path.abspath(self.infile))
        filename = os.path.basename(os.path.abspath(self.infile))
        cmd = f"python {scfea} --data_dir {self.refdir} --moduleGene_file {self.moduleGene} --stoichiometry_matrix {self.stoichiometry} --sc_imputation {self.imputation} --input_dir {indir} --test_file {filename} --res_dir {self.outdir} "
        if self.cName:
            cmd += f" --cName_file {self.cName}"
        logging.info(cmd)
        print(cmd)
        subprocess.run(cmd, shell=True)

def run_scfea_cli():
    """cell-wise metabolic estimation using scFEA in one command

    Online Manual:
    - https://github.com/changwn/scFEA/blob/master/scFEA_tutorial1.ipynb
    - https://github.com/changwn/scFEA/blob/master/scFEA_tutorial2.ipynb

    Requirements:
    - scFEA.py: add to PATH, https://github.com/changwn/scFEA

    Attributes:
        infile (str): Input seurat count matrix csv file
        outdir (str): Output directory
        refdir (str): directory including scFEA model files
        species (str): human or mouse
        moduleGene (str): table contains genes for each module, default is module_gene_m168.csv for human
        stoichiometry(str): table contains relations between compounds and modules, default is cmMat_171.csv for human
        cName(str): The name of compounds, 2 rows: compounds name and id
        imputation(bool): Whether perform imputation for SC dataset

    Raises:
        ModuleNotFoundError Need install scFEA first!!!
        OSError Need different files to run!!!
    """
    import argparse
    parser = argparse.ArgumentParser( description='cell-wise metabolic estimation using scFEA ')
    parser.add_argument( 'infile', type=str, help='Input seurat count matrix csv file')
    parser.add_argument( 'outdir', type=str, help='output directory')
    parser.add_argument( '--refdir', type=str, help='scFEA model directory')
    parser.add_argument( '--species', type=str, default = "human",
                        help='human or mouse, default is human')
    parser.add_argument( '--moduleGene', type=str, help='table contains genes for each module')
    parser.add_argument( '--stoichiometry', type=str, help='table contains relations between compounds and modules')
    parser.add_argument( '--cName', type=str, help='The name of compounds, 2 rows: compounds name and id')
    parser.add_argument( '--imputation', type=bool, default = True,
                        help='Whether perform imputation for SC data, default is True')
    args = parser.parse_args()
    args = vars(args)
    print(args)
    obj = MetabolicFluxUsingSCFEA(**args)
    obj.run()

from setuptools import setup, find_packages
setup(
    name="pyscrnarecom",
    version="0.1.0",
    author="SuMin",
    author_email="sumin2012@163.com",
    description=" ...",
    url="https://github.com/WXlab-NJMU/scrna-recom/",
    python_requires=">=3.10",
    license = "MIT",

    # package data specified in MANIFEST.in
    include_package_data=True,
    exclude_package_data={'':['.gitignore']},
    packages=['pyscrnarecom'],
    # package structure: source path, target path
    package_dir = {"pyscrnarecom": "src/pyscrnarecom"},

    # package will executable commands
    #   'cli-name = mypkg.mymodule:some_func',
    entry_points={
        'console_scripts': [
            'scrna-basic=pyscrnarecom.processing:run_scanpy_cli',
            'scrna-basic-plot-focus=pyscrnarecom.processing:plot_genes_cli',
            'scrna-rawdata=pyscrnarecom.rawdata:run_cellranger_cli',
            'scrna-trajectory=pyscrnarecom.trajectory:run_scvelo',
            'scrna-regulon=pyscrnarecom.network:run_pyscenic_cli',
            'scrna-metabolism=pyscrnarecom.metabolism:run_scfea_cli'
        ]
    },

    # required packages
    ## platfrom specific dependencies:
    ##   "pywin32 >= 1.0;platform_system=='Windows'"
    install_requires=[
        "Cython",
        "pytz",
        "numpy",
        "docutils",
        "sphinx",
        "sphinx-book-theme",
        "scvelo",
        "velocyto",
        "pyscenic",
        "scanpy",
        "loompy",
        "matplotlib",
        "seaborn",
        "multicoretsne",
        "umap-learn",
        "anndata",
    ],
    # extra requirements: identifier: [required_packages]
    extras_require={
        "PDF": ["ReportLab>=1.2", "RXP"],
    },
)

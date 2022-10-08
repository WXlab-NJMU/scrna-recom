from setuptools import setup, find_packages
setup(
    name="pyscrnarecom",
    version="1.0",
    author="SuMin",
    author_email="sumin2012@163.com",
    description=" ...",
    url="http://github.com/",
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
            'basic-scanpy=pyscrnarecom.processing:run_scanpy_cli',
            'basic-plot-focus=pyscrnarecom.processing:plot_genes_cli',
            'basic-rawdata=pyscrnarecom.rawdata:run_cellranger_cli',
            'trajectory-scvelo=pyscrnarecom.trajectory:run_scvelo',
            'network-pysenic=pyscrnarecom.network:run_pyscenic_cli',
            'metabolism-scfea=pyscrnarecom.metabolism:run_scfea_cli'
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
        "MulticoreTSNE",
        "umap-learn",
    ],
    # extra requirements: identifier: [required_packages]
    extras_require={
        "PDF": ["ReportLab>=1.2", "RXP"],
    },
)

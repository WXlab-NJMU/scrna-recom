from setuptools import setup, find_packages
setup(
    name="scrna",
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
    # package structure: source path, target path
    package_dir = {
        ".": "src"
    },

    # package will executable commands
    #   'cli-name = mypkg.mymodule:some_func',
    entry_points={
        'console_scripts': [
            'run-scvelo=scrna.trajectory:run_scvelo',
            'run-scanpy=scrna.processing:run_scanpy_cli',
            'plot-scanpy-focus=scrna.processing:plot_genes_cli',
            'run-pysenic=scrna.network:run_pyscenic_cli',
        ]
    },

    # required packages
    ## platfrom specific dependencies:
    ##   "pywin32 >= 1.0;platform_system=='Windows'"
    install_requires=[
        "docutils",
        "sphinx",
        "sphinx-book-theme",
        "scvelo",
        "velocyto",
        "pyscenic",
        "scanpy",
        "numpy",
        "loompy",
        "matplotlib",
        "seaborn",
    ],
    # extra requirements: identifier: [required_packages]
    extras_require={
        "PDF": ["ReportLab>=1.2", "RXP"],
    },
)

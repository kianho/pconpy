from setuptools import setup

setup(name="pconpy",
    version="0.1",
    description="a python module for generating protein contact maps and distance maps",
    author="Kian Ho",
    author_email="hui.kian.ho@gmail.com",
    url="http://www.github.com/kianho/pconpy",
    license="MIT",
    packages=["pconpy"],
    keywords=["protein", "protein structure", "bioinformatics", "contact map",
        "distance map", "computational biology", "visualization",
        "visualisation", "amino acids"],

    # Dependencies
    #install_requires=["biopython>=1.65", "numpy>=1.9.0", "matplotlib>=1.4.0", "docopt>=0.6.2"]
    #install_requires=["biopython>=1.65", "matplotlib>=1.4.0", "docopt>=0.6.2"]
    install_requires=["docopt>=0.6.2"]
)

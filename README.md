PConPy—a Python module for generating 2D protein maps
=====================================================

This is the official repository for the redevelopment of PConPy, previously
published as:

    Hui Kian Ho, Michael J. Kuiper, and Kotagiri Ramamohanarao (2008). “PConPy–a
    Python module for generating 2D protein maps”. In: Bioinformatics 24,
    pp. 2934-2935.

The full article is accessible
[here](http://bioinformatics.oxfordjournals.org/content/24/24/2934.full).

The original source code associated with the article can be accessed from the
`legacy` branch of this repository.


## Motivation

_Contact maps_ are a reduced representation of 3D structure that are invariant
to affine transformations. This enables the researcher to perform a rapid
visual exploratory data analysis of a protein structure without the need for 3D
rendering. Additionally, the 2D matrix representation of a contact map (known
as a _contact matrix_) can naturally be used as numeric input in subsequent
automated knowledge discovery or machine learning tasks.

## History

This project was conceived at the beginning of my PhD when I first became
interested in analysing [known 3D protein structures](http://www.pdb.org) for
interesting patterns. The project also provided me with one of my first
opportunities to experiment with the [matplotlib](http://matplotlib.org/) library.


## Acknowledgements

Peter Cock's (@peterjc) [tutorial](http://goo.gl/q7DNt7) was especially helpful
in getting me started with the mechanics of PDB file parsing and visualisation.

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
visual exploratory data analysis of a protein structure without the need for
rendering or geometric operations in 3D. Additionally, the 2D matrix
representation of a contact map (known as a _contact matrix_) can be easily
manipulated for use by conventional machine learning algorithms for protein
structure prediction tasks.

## History

This project was conceived at the beginning of my PhD when I first became interested in
analysing [known 3D protein structures](http://www.pdb.org) for interesting
patterns to be (potentially) used for protein structure prediction tasks.
This project provided me with one of my first opportunities to experiment
with the [matplotlib]() library.


## Acknowledgements

Peter Cock's (@peterjc) [tutorial](http://goo.gl/q7DNt7) was especially helpful
in getting me started with the mechanics of PDB file parsing and visualisation.

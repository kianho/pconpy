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


## About

A [_contact map_](http://en.wikipedia.org/wiki/Protein_contact_map) is a 2D
representation of protein structure that illustrates the presence or absence of
contacts between individual amino acids. This enables the rapid visual
exploratory data analysis of structural features without 3D rendering software.
Additionally, the underlying 2D matrix of a contact map, known as a _contact
matrix_, can naturally be used as numeric input in subsequent automated
knowledge discovery or machine learning tasks. PConPy generates
publication-quality renderings of contact maps, and the related distance- and
hydrogen bond maps.

## History

This project was conceived at the beginning of my PhD when I first became
interested in analysing [known 3D protein structures](http://www.pdb.org) for
interesting patterns. The project also provided me with one of my first
opportunities to experiment with the [matplotlib](http://matplotlib.org/) library.


## Acknowledgements

Peter Cock's (@peterjc) [tutorial](http://goo.gl/q7DNt7) was especially helpful
in getting me started with the mechanics of PDB file parsing and visualisation.

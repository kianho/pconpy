#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    pconpy.py -p PDB [-c CHAINID] -o OUTPUT [options]

Options:
    -p PDB, --pdb PDB
    -c CHAINIDS    Space-separated list of chain identifiers (defaults to the first chain if left empty).
    -o OUTPUT, --output OUTPUT
    -m DIST                     The inter-residue distance measure (see below) [default: CA].
    -t THRESH                   The contact distance threshold.
    --plaintext                 Generate a plaintext distance/contact matrix
                                and write to stdout (recommended for
                                piping into other CLI programs).
    --title TITLE               The title of the plot (optional).
    --font FONT                 Set the font family (via matplotlib).
    --noframe
    -v, --verbose

Distance measures (-m DIST):
    "ca" -- Conventional CA-CA distance, this is the default distance measure.
    "cb" -- The CB-CB distance.
    "cmass" -- The distance between the residue centers of mass.
    "sccmass" -- The distance between the sidechain centers of mass
    "minvdw" -- The minimum distance between the VDW radii of each residue.

Dependencies:
    docopt
    biopython
    numpy
    matplotlib/pylab

To do:
    - centres-of-mass distance
    - optional colour bar for distance map
    - greyscale distance map

"""

import mplutils
import matplotlib as mpl
import os
import sys
import re
import numpy
import pylab
import tempfile

from collections import defaultdict
from itertools import ifilter, product
from itertools import combinations, combinations_with_replacement
from docopt import docopt

import Bio.PDB

PWD = os.path.dirname(os.path.abspath(__file__))

# The atom names of the backbone and sidechain atoms are based on those defined
# in the `ProDy` module.
BACKBONE_ATOMS = set(["CA", "C", "O", "N"])
BACKBONE_FULL_ATOMS = set(['CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3', 'OXT'])

# VDW radii values
#
# NOTE: 
#   - The first character in the PDB atom name is the element of the atom.
#
VDW_RADII = { "N" : 1.55, "C" : 1.70,"H" : 1.20, "O" : 1.52, "S": 1.85 }


def is_backbone(atom):
    return atom.get_id() in BACKBONE_ATOMS

def is_sidechain(atom):
    return atom.get_id() not in BACKBONE_FULL_ATOMS

def get_backbone_atoms(res):
    return filter(lambda atom : is_backbone(atom), res.get_iterator())

def get_sidechain_atoms(res, infer_CB=True):
    """

    Arguments:
        res --
        infer_CB -- infer a virtual CB atom if the residue does not contain
        sidechain atoms e.g. glycine residues (default: True).

    Returns:
        ...

    """
    atoms = filter(lambda atom : is_sidechain(atom), res.get_iterator())

    if (not atoms) and infer_CB:
        CB_coord = get_atom_coord(res, "CB")
        CB_atom = Bio.PDB.Atom.Atom("CB", CB_coord,
                    None, None, None, None, None, "C")
        atoms = [CB_atom]

    return atoms


def pdb_to_ag(pdb_file, chain_id=None):
    """Generate the AtomGroup object from a single PDB file.

    """

    dssp_fd, dssp_out = tempfile.mkstemp()
    dssp = pd.proteins.execDSSP(pdb_file, outputname=dssp_out)

    ag = pd.proteins.parsePDB(pdb_file)
    ag = pd.proteins.parseDSSP(dssp, ag, parseall=True)

    # chain = ag.select("protein").getHierView()[opts["--chid"]]

    os.close(dssp_fd)
    os.remove(dssp_out)

    return ag.select("protein")


def get_residues(pdb_fn, chain_ids=None, model_num=0):
    """Build a simple list of residues from a single chain.

    Arguments:
        pdb_fn --
        chain_ids -- (default: None).
        model_num -- (default: 0).

    Returns:
        ...

    """

    pdb_id = os.path.splitext(os.path.basename(pdb_fn))[0]

    parser = Bio.PDB.PDBParser(pdb_id, pdb_fn)
    struct = parser.get_structure(pdb_id, pdb_fn)
    model = struct[model_num]

    if chain_ids is None:
        # get residues from every chain.
        chains = model.get_list()
    else:
        chains = [ model[ch_id] for ch_id in chain_ids ]

    residues = []

    # To do, handle:
    # - disordered residues/atoms (todo)
    # - non-amino acids
    # - non-standard amino acids
    for ch in chains:
        for res in ifilter(lambda r : Bio.PDB.is_aa(r), ch.get_residues()):
            if not Bio.PDB.is_aa(res, standard=True):
                sys.stderr.write("WARNING: non-standard AA at %r%s"
                        % (res.get_id(), os.linesep))
            residues.append(res)

    return residues


def calc_eucl_distance(res_a, res_b, atom_a="CA", atom_b="CA"):
    """Compute the Euclidean distance between specific atoms of a pair of
    residues.

    Arguments:
        res_a --
        res_b --
        atom_a --
        atom_b --

    Returns:
        ...

    """

    A = res_a[atom_a].get_coord()
    B = res_b[atom_b].get_coord()

    return numpy.linalg.norm(A, B)


def get_atom_coord(res, atom_name, verbose=False):
    """
    If the CB coordinate exists, return it, otherwise compute a
    virtual CB coordinate by rotating the N atom -120 degrees around the CA-C
    vector. This will occur with Glycine residues which don't have sidechain
    atoms.

    Arguments:
        res --
        atom_name --
        verbose --

    Reference:
        http://goo.gl/OaNjxe

    """

    try:
        coord = res[atom_name].get_coord()
    except KeyError:
        if atom_name != "CB":
            raise NotImplementedError

        if verbose:
            sys.stderr.write(
                "WARNING: computing virtual {} coordinate.".format(atom_name)
                    + os.linesep)

        assert("N" in res)
        assert("CA" in res)
        assert("C" in res)

        # NOTE:
        # These are Bio.PDB.Vector objects and _not_ numpy arrays.
        N = res["N"].get_vector()
        CA = res["CA"].get_vector()
        C = res["C"].get_vector()

        CA_N = N - CA
        CA_C = C - CA

        rot_mat = Bio.PDB.rotaxis(numpy.pi * 120.0 / 180.0, CA_C)

        coord = (CA + N.left_multiply(rot_mat)).get_array()

    return coord


def calc_center_of_mass(atoms):
    """Compute the center of mass from a list of atoms.

    Arguments:
        atoms -- a list of Bio.PDB.Atom.Atom objects.

    Returns:
        the center of mass.

    """
    return ( numpy.array([(x.mass * x.get_coord()) for x in atoms]).sum()
                / sum(x.mass for x in atoms) )


def calc_minvdw_distance(res_a, res_b):
    """

    """

    min_dist = None

    for a in res_a.get_iterator():
        for b in res_b.get_iterator():
            radii_a = VDW_RADII.get(a.get_id()[0], 0.0)
            radii_b = VDW_RADII.get(b.get_id()[0], 0.0)

            A = a.get_coord()
            B = b.get_coord()

            dist = numpy.linalg.norm(A - B) - radii_a - radii_b

            if (min_dist is None) or dist < min_dist:
                min_dist = dist

    return min_dist


def calc_cmass_distance(res_a, res_b, sidechain_only=False):
    """Compute the distance between the centres of mass of both residues.

    Arguments:
        res_a --
        res_b --
        sidechain_only --

    Returns:
        ...

    """

    if sidechain_only:
        atoms_a = get_sidechain_atoms(res_a, infer_CB=True)
        atoms_b = get_sidechain_atoms(res_b, infer_CB=True)
    else:
        atoms_a = res_a.get_list()
        atoms_b = res_b.get_list()

    A = calc_center_of_mass(atoms_a)
    B = calc_center_of_mass(atoms_b)

    return numpy.linalg.norm(A-B)


def calc_distance(res_a, res_b, metric="CA"):
    """Calculate the Euclidean distance between a pair of residues according to
    a given distance metric.

    Arguments:
        res_a --
        res_b --
        metric -- the distance metric (default: "CA").

    Returns:
        ...

    """

    if metric in ("CA", "CB"):
        A = get_atom_coord(res_a, metric)
        B = get_atom_coord(res_b, metric)

        dist = numpy.linalg.norm(A-B)
    elif metric == "cmass":
        dist = calc_cmass_distance(res_a, res_b)
    elif metric == "sccmass":
        dist = calc_cmass_distance(res_a, res_b, sidechain_only=True)
    elif metric == "minvdw":
        dist = calc_minvdw_distance(res_a, res_b)
    else:
        raise NotImplementedError

    return dist


def calc_dist_matrix(residues, metric="CA", threshold=None, symmetric=False):
    """Calculate the distance matrix for a list of residues.

    Arguments:
        residues --
        metric --
        threshold -- (default: None)
        symmetric -- (default: False)

    Returns
        ...

    """

    mat = numpy.zeros((len(residues), len(residues)), dtype="float64")

    # after the distances are added to the upper-triangle, the nan values
    # indicate the lower matrix values, which are "empty", but can be used to
    # convey other information if needed.
    mat[:] = numpy.nan

    # Compute the upper-triangle of the underlying distance matrix.
    #
    # TODO:
    # - parallelise this over multiple processes + show benchmark results.
    # - use the lower-triangle to convey other information.
    pair_indices = combinations_with_replacement(xrange(len(residues)), 2)

    for i, j in pair_indices:
        res_a = residues[i]
        res_b = residues[j]
        dist = calc_distance(res_a, res_b, metric)

        mat[i,j] = dist

        symmetric = True
        if symmetric:
            mat[j,i] = dist

    # we transpose i with j so the distances are contained only in the
    # upper-triangle.
    mat = mat.T
    mat = numpy.ma.masked_array(mat, numpy.isnan(mat))

    if threshold:
        mat = numpy.ma.masked_greater_equal(mat, threshold)

    return mat


def mat_to_ascii(mat):
    """
    """

    if "float" in mat.dtype.name:
        fmt = lambda x : "%.2f" % x
    else:
        fmt = lambda x : ( "%d" % x if x > 0 else " " )

    for row in mat:
        print "".join( fmt(val) for val in row )

    return


if __name__ == '__main__':
    opts = docopt(__doc__)

    if opts["-t"]:
        opts["-t"] = float(opts["-t"])

    if opts["-c"]:
        chain_ids = opts["-c"].upper().split(",")

        # Check that pdb chain ids are alphanumeric (see:
        # http://deposit.rcsb.org/adit/).
        if not numpy.all(map(str.isalnum, chain_ids)):
            sys.stderr.write()

    residues = get_residues(opts["--pdb"])

    #
    # Generate the underlying 2D matrix for the selected plot.
    #
    mat = calc_dist_matrix(residues, metric=opts["-m"],
            threshold=opts["-t"])

    if opts["--plaintext"]:
        if opts["-t"] is not None:
            # Use mask distances below the selected threshold.
            mat = mat < float(opts["-t"])
            fmt = "%d"
        else:
            fmt = "%.3f"

        numpy.savetxt(sys.stdout, mat, fmt=fmt)
    else:
        mplutils.init_pylab(family=opts["--font"])

        # hide all the spines i.e. no axes are drawn
        mplutils.init_spines(hidden=["top", "bottom", "left", "right"])

        # make figure square-shaped
        pylab.gcf().set_figwidth(6.0)
        pylab.gcf().set_figheight(6.0)

        pylab.xlim(0, len(residues))
        pylab.ylim(0, len(residues))

        pylab.xlabel("Residue index")
        pylab.ylabel("Residue index")

        ax, fig = pylab.gca(), pylab.gcf()

        if opts["-t"] is None:
            map_obj = pylab.pcolormesh(mat, shading="flat", edgecolors="None", cmap=mpl.cm.jet_r)

            # draw the colour bar
            box = ax.get_position()
            pad, width = 0.02, 0.02
            cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
            cbar = pylab.colorbar(map_obj, drawedges=False, cax=cax)
            cbar.outline.set_visible(False)
            pylab.ylabel("Distance (Angstroms)")
        else:
            if not opts["--noframe"]:
                mplutils.init_spines(hidden=[])
            map_obj = pylab.pcolormesh(mat.astype("int"),
                shading="flat", edgecolors="None", cmap=mpl.cm.Greys)

        if opts["--title"] is not None:
            ax.set_title(opts["--title"], fontweight="bold")

        pylab.savefig(opts["--output"], bbox_inches="tight")

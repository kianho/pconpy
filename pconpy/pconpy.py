#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    pcon.py -p PDB [-c CHAINID] -o OUTPUT [options]

Options:
    -p PDB, --pdb PDB
    -c CHAINIDS, --chid CHAINIDS    Space-separated list of chain identifiers (defaults to the first chain if left empty).
    -o OUTPUT, --output OUTPUT
    -d DIST                     The inter-residue distance measure (see below) [default: ca].
    -t THRESH                   The contact distance threshold [default: 8.0].
    --plaintext                 Generate a plaintext distance/contact matrix
                                and write to stdout (recommended for
                                piping into other CLI programs).
    --title TITLE               The title of the plot (optional).
    --font FONT                 Set the font family (via matplotlib).
    --noframe
    -v, --verbose

Distance measures (-d DIST):
    "ca" -- Conventional CA-CA distance, this is the default distance measure.
    "cb" -- The CB-CB distance.
    "cmass" -- The distance between the residue centers of mass.
    "sccmass" -- The distance between the sidechain centers of mass
    "minvdw" -- The minimum distance between the VDW radii of each residue.

Dependencies:
    docopt
    prody
    numpy
    matplotlib/pylab

To do:
    - centres-of-mass distance
    - optional colour bar for distance map
    - greyscale distance map

"""

import mplutils
import matplotlib
import matplotlib as mpl
import os
import sys
import numpy
import numpy as np
import pylab
import tempfile
import Bio.PDB

from docopt import docopt
from collections import defaultdict
from itertools import ifilter, combinations

VDW_RADII = \
    defaultdict(float,
        { "N" : 1.55, "C" : 1.70,"H" : 1.20, "O" : 1.52, "S": 1.85 })


import Bio.PDB


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


def get_residues(pdb_fn, chain_id=None, model_num=0):
    """Build a simple list of residues from a single chain.

    Arguments:
        pdb_fn --
        chain_id --
        model_num --

    Returns:
        ...

    """

    pdb_id = os.path.splitext(os.path.basename(pdb_fn))[0]

    parser = Bio.PDB.PDBParser(pdb_id, pdb_fn)
    struct = parser.get_structure(pdb_id, pdb_fn)
    model = struct[model_num]

    if chain_id is None:
        # get residues from every chain.
        chains = model.get_list()
    else:
        chains = [model[chain_id]]

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


def calc_distance(res_a, res_b, metric="ca"):
    """Calculate the Euclidean distance between a pair of residues according to
    a given distance metric.

    """

    if metric in ("ca", "cb"):
        # Conventional CA-CA or CB-CB distance.
        atom_name = metric.upper()

        assert(atom_name in res_a)
        assert(atom_name in res_b)

        atom_a = res_a[atom_name] 
        atom_b = res_b[atom_name]

        A = atom_a.get_coord()
        B = atom_b.get_coord()

        dist = numpy.linalg.norm(A-B)
    else:
        raise NotImplementedError

    return dist


def make_dist_matrix(residues, metric="ca"):
    """Calculate the distance matrix for a list of residues.

    Arguments:
        residues --
        metric --

    Returns
        ...

    """

    mat = numpy.zeros((len(residues), len(residues)), dtype="float64")

    # Compute the upper-triangle of the underlying distance matrix.
    #
    # TODO:
    # - parallelise this over multiple processes + show benchmark results.
    for i, j in combinations(xrange(len(residues)), 2):
        res_a = residues[i]
        res_b = residues[j]

        mat[i,j] = calc_distance(res_a, res_b, metric)

    return mat



def calc_minvdw_distance(res_a, res_b):
    """
    """

    atoms_a = list(res_a.iterAtoms())
    atoms_b = list(res_b.iterAtoms())
    min_dist = None

    for a in atoms_a:
        for b in atoms_b:
            dist = pd.measure.calcDistance(a, b)

            if (min_dist is None) or dist < min_dist:
                min_dist = dist

    return min_dist


def calc_cmass_distance(res_a, res_b, sidechain_only=False):
    """
    """
    _res_a = res_a
    _res_b = res_b

    if sidechain_only:
        res_a = res_a.select("sidechain")
        res_b = res_b.select("sidechain")

    if res_a is None:
        res_a = _res_a

    if res_b is None:
        res_b = _res_b

    coord_a = pd.measure.calcCenter(res_a, weights=res_a.getMasses())
    coord_b = pd.measure.calcCenter(res_b, weights=res_b.getMasses())

    return pd.measure.calcDistance(coord_a, coord_b)


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

    residues = get_residues(opts["--pdb"])

    #
    # Generate the underlying 2D matrix for the selected plot.
    #

    if opts["-d"] in ("ca", "cb"):
        # distance between specific atoms.
        atom_name = opts["-d"].upper()
        res_coords = []

        for res in residues:
            try:
                coord = res[atom_name].getCoords()
            except KeyError:
                coord = res["CA"].getCoords()
            res_coords.append(coord)

        res_coords = np.array(res_coords)
        mat = pd.measure.buildDistMatrix(res_coords, res_coords)
    else:
        mat = np.zeros((len(residues), len(residues)), dtype="float64")

        if opts["-d"] == "minvdw" :
            for i, res_a in enumerate(residues):
                for j, res_b in enumerate(residues):
                    mat[i,j] = calc_minvdw_distance(res_a, res_b)
        elif opts["-d"] == "cmass":
            for i, res_a in enumerate(residues):
                for j, res_b in enumerate(residues):
                    mat[i,j] = calc_cmass_distance(res_a, res_b)
        elif opts["-d"] == "sccmass":
            empty_indices = []
            for i, res_a in enumerate(residues):
                for j, res_b in enumerate(residues):
                    dist = calc_cmass_distance(res_a, res_b, sidechain_only=True)
                    if dist is None:
                        empty_indices.append((i,j))
                        continue
                    mat[i,j] = dist

            for i, j in empty_indices:
                mat[i,j] = mat.max()

    if opts["--plaintext"]:
        if opts["-t"] is not None:
            mat = mat < float(opts["-t"])
            fmt = "%d"
        else:
            fmt = "%.3f"

        np.savetxt(sys.stdout, mat, fmt=fmt)
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
            thresh = float(opts["-t"])
            if not opts["--noframe"]:
                mplutils.init_spines(hidden=[])
            map_obj = pylab.pcolormesh((mat < thresh).astype("int"),
                shading="flat", edgecolors="None", cmap=mpl.cm.Greys)

        if opts["--title"] is not None:
            ax.set_title(opts["--title"], fontweight="bold")

        pylab.savefig(opts["--output"], bbox_inches="tight")

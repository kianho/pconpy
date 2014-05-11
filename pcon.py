#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    pcon.py -p PDB [-c CHAINID] [-o OUTPUT] [-d DIST] [-t THRESH]

Options:
    -p PDB, --pdb PDB
    -c CHAINID, --chid CHAINID
    -o OUTPUT, --output OUTPUT
    -d DIST                         [default: CA]
    -t THRESH
    -v, --verbose

Dependencies:
    docopt
    prody
    numpy
    matplotlib/pylab

To do:
    - centres-of-mass distance
    - min. VDW distance
    - optional colour bar for distance map
    - greyscale distance map

"""

import mplutils
import matplotlib as mpl
import os
import sys
import numpy as np
import prody as pd
import pylab
import tempfile

from docopt import docopt
from collections import defaultdict

VDW_RADII = \
    defaultdict(float,
        { "N" : 1.55, "C" : 1.70,"H" : 1.20, "O" : 1.52, "S": 1.85 })


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


def get_residues(pdb_file, chain_id=None):
    """
    """

    ag = pdb_to_ag(pdb_file)

    # Default to the first/only chain ID if none were specified.
    if chain_id is None:
        chain_id = ag.getChids()[0]

    chain = ag.getHierView()[chain_id]

    return chain.iterResidues()


def calc_minvdw_distance(res_a, res_b):
    """
    """

    atoms_a = list(res_a.iterAtoms())
    atoms_b = list(res_b.iterAtoms())
    min_dist = None

    for a in atoms_a:
        for b in atoms_b:
            dist = pd.measure.calcDistance(a, b)

            if ( min_dist is None ) or dist < min_dist:
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

    # hide all the spines i.e. no axes are drawn
    mplutils.init_spines(hidden=["top", "bottom", "left", "right"])

    # make figure square-shaped
    pylab.gcf().set_figwidth(6.0)
    pylab.gcf().set_figheight(6.0)

    residues = list(get_residues(opts["--pdb"]))

    if opts["-d"] in [ "CA" ]:
        res_coords = np.array([ r[opts["-d"]].getCoords() for r in residues ])
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
        elif opts["-d"] == "scmass":
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

    pylab.xlim(0, len(residues))
    pylab.ylim(0, len(residues))

    if opts["-t"] is None:
        map_obj = pylab.pcolormesh(mat, shading="flat", edgecolors="None", cmap=mpl.cm.jet_r)

        # draw the colour bar
        ax, fig = pylab.gca(), pylab.gcf()
        box = ax.get_position()
        pad, width = 0.02, 0.02
        cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
        cbar = pylab.colorbar(map_obj, drawedges=False, cax=cax)
        cbar.outline.set_visible(False)
    else:
        thresh = float(opts["-t"])
        mplutils.init_spines(hidden=[])
        map_obj = pylab.pcolormesh((mat < thresh).astype("int"),
            shading="flat", edgecolors="None", cmap=mpl.cm.Greys)

    pylab.savefig(opts["--output"], bbox_inches="tight")

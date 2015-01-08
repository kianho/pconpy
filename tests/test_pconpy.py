#!/usr/bin/env python
# encoding: utf-8
"""

Date:
    Wed Jan  7 23:58:53 AEDT 2015

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    test_pconpy.py

"""

import os
import sys

#import pconpy
import numpy
import unittest
import Bio.PDB

from itertools import ifilter, combinations

PWD = os.path.dirname(os.path.realpath(__file__))
PDB_DIR = os.path.join(PWD, "pdb_files")


def load_pdb(pdb_fn, chain_id=None):
    """Load a single chain from a pdb file.

    """

    return pconpy.pdb_to_ag(pdb_fn, chain_id)


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
        chain = model.get_list()[0]
    else:
        chain = model[chain_id]

    residues = []

    # To do, handle:
    # - disordered residues/atoms (todo)
    # - non-amino acids
    for res in ifilter(lambda r : Bio.PDB.is_aa(r), chain.get_residues()):
        if not Bio.PDB.is_aa(res):
            sys.stderr.write("WARNING: non-standard AA at %r%s"
                    % (res.get_id(), os.linesep))
        residues.append(res)

    return residues


def calc_distance(res_a, res_b, method="ca"):
    """Calculate the Euclidean distance between a pair of residues according to
    a specified distance method.

    """

    if method in ("ca", "cb"):
        # Conventional CA-CA or CB-CB distance.
        atom_name = method.upper()

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


class TestGeometry(unittest.TestCase):
    def setUp(self):
        """

        """

        pdb_id = "1ejg"
        pdb_fn = os.path.join(PDB_DIR, "{}.pdb".format(pdb_id))

        self.residues = get_residues(pdb_fn, "A")

        return


    def tearDown(self):
        """

        """

        return


    def test_calc_distance(self):
        """

        """

        #pconpy.calc_eucl_distance()
        residues = self.residues

        # upper triangle indices only
        for i, j in combinations(range(len(residues)), 2):
            res_i = residues[i]
            res_j = residues[j]

            print calc_distance(res_i, res_j, "ca")


        return

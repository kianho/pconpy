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

import pconpy
import prody
import unittest

PWD = os.path.dirname(os.path.realpath(__file__))
PDB_DIR = os.path.join(PWD, "pdb_files")

ATOM_NAMES


def load_pdb(pdb_fn, chain_id=None):
    """Load a single chain from a pdb file.

    """

    return pconpy.pdb_to_ag(pdb_fn, chain_id)


class TestGeometry(unittest.TestCase):
    def setUp(self):
        """

        """

        pdb_fn = os.path.join(PDB_DIR, "1ubq.pdb")

        self._residues = list(pconpy.get_residues(pdb_fn, "A"))

        return


    def tearDown(self):
        """

        """

        return


    def test_calc_eucl_distance(self):
        """

        """

        pconpy.calc_eucl_distance()

        return

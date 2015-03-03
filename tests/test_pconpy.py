#!/usr/bin/env python
# encoding: utf-8
"""

Date:
    Wed Jan  7 23:58:53 AEDT 2015

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    This file contains the unit tests for pconpy.

Usage:
    This script can be run by calling "nosetests" or "py.test" from the root
    pconpy directory.

"""

import os
import sys

import pconpy
import numpy

PWD = os.path.dirname(os.path.realpath(__file__))
PDB_DIR = os.path.join(PWD, "pdb_files")


class TestGeometry:
    # Skeleton fixtures (will use in the future)
    @classmethod
    def setup_class(self):
        return

    @classmethod
    def teardown_class(self):
        return

    def pdb_id_to_fn(self, pdb_id):
        """Get the path to a test pdb file using its pdb id.

        """
        return os.path.join(PDB_DIR, "{}.pdb".format(pdb_id))


    def test_calc_dist_matrix(self):
        """Test the generation of distance matrices.

        """

        def run_test(metric):
            """
            """

            # Ideal case should result in a square non-zero distance matrix.
            pdb_fn = self.pdb_id_to_fn("1ubq")
            residues = pconpy.get_residues(pdb_fn)
            mat = pconpy.calc_dist_matrix(residues, metric)

            assert( mat is not None )
            assert( numpy.any(mat > 0) )
            assert( mat.shape == (len(residues), len(residues)) )

            # Empty residues list should result in an empty distance matrix.
            residues = []
            mat = pconpy.calc_dist_matrix(residues, metric=metric)

            assert( mat is not None )
            assert( mat.shape == (0,0) )

            return

        # Test each of the distance metrics on the ideal and edge cases.
        for metric in ("CA", "CB", "cmass", "sccmass", "minvdw"):
            yield run_test, metric

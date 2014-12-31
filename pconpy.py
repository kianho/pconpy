#!/usr/bin/env python
"""pconpy.py - generate 2D protein maps from PDB files.

File: pconpy.py
Author: Hui Kian Ho
Created: July 2008

Description:
    PConPy is an open-source Python module for generating 2D protein maps,
    these include contact maps, distance maps, and hydrogen-bond plots.

Dependencies:
    PConPy depends on a number of Python modules which are all available via
    the Synaptic Package Manager on Ubuntu Linux:
        - NumPy (http://www.numpy.org)
        - ScientificPython
            (http://dirac.cnrs-orleans.fr/plone/software/scientificpython/)
            Note that this is *not* the same as SciPy.
        - Matplotlib (http://matplotlib.sourceforge.net/)
            The pylab module is included in Matplotlib.

    PConPy depends on the following external programs:
        - STRIDE (http://webclu.bio.wzw.tum.de/stride/)
        - DSSP (http://swift.cmbi.ru.nl/gv/dssp/)

Comments:
    I've kept all classes and functions inside a single file to simplify
    portability.

    More plotting functionality will be added as this project progresses.

Contributors:
    KH - Kian Ho (author of pconpy)
    MK - Matthew Kallada (pconpy user)

"""

import os
import sys
import tempfile
import traceback


# Check for external Python modules.
try:
    import pylab
except ImportError:
    sys.stderr.write(
        "Error: Matplotlib is not installed, it can be downloaded " + \
        "from http://matplotlib.sourceforge.net/\n"
    )

try:
    import Scientific
except ImportError:
    sys.stderr.write(
        "error: ScientificPython is not installed, it can be downloaded " + \
        "from http://dirac.cnrs-orleans.fr/plone/software/scientificpython/\n"
    )

try:
    import numpy
except ImportError:
    sys.stderr.write(
        "error: NumPy is not installed, it can be downloaded " + \
        "from http://www.numpy.org/\n"
    )


from Scientific.IO.PDB import Structure
from optparse import OptionParser
from numpy import arccos, arange, array, cross, dot
from numpy import matrix, pi, sin, sqrt, zeros
from numpy.linalg import norm


# Specify the default temp file directory.
tempfile.tempdir = "/var/tmp/"

# Non-side-chain atoms.
representative_atoms = ["CA", "CB", "N", "C", "O"]

# Lookup table to convert 3-letter AA names to 1-letter names.
aa_index_3_to_1 = { "ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
                    "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
                    "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
                    "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
                    "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y",
                    "UNK" : "X" }

# Atomic masses used to calculate the center of mass of each residue.
# Taken from http://en.wikipedia.org/wiki/Atomic_mass
MASS = { "N" : 14.01, "C" : 12.01, "H" : 1.008, "O" : 16.00, "S" : 32.07 }

# VDW radii of atoms in Angstroms.
# Taken from http://www.imb-jena.de/ImgLibDoc/glossary/IMAGE_VDWR.html
VDW_RADII = { "N" : 1.55, "C" : 1.70,"H" : 1.20, "O" : 1.52, "S": 1.85 }

# SSE colouring scheme.
SSE_COLOUR = {"H" : "m",        # alpha-helix --> magenta
              "E" : "#ffd700",  # beta-strand --> yellow
              "B" : "orange",   # isolated beta-bridge --> orange
              "b" : "orange",   # isolated beta-bridge --> orange
              "T" : "#30D5C8",  # turn --> cyan
              "C" : "#D2B48C",  # coil --> light brown
              "G" : "b",        # 3/10-helix --> blue
              "I" : "g",        # pi-helix --> green
              "S" : "#d2b48c",  # bend --> light brown
              None : "black"}

# Name of the font to use for the axes labels and plot title.
FONT_NAME = "Arial"

FAKE_CHID = "?"


def default_chain_id(chain_id, default=FAKE_CHID):
    """Return a "default" chain identifier if the raw parsed chain id (from dssp
    or stride) is not an alpha-numeric character.

    """

    return chain_id if chain_id.isalnum() else default


class Protein(object):
    """A wrapper class around Scientific.IO.PDB.Structure to simplify access
    to protein residues and chains.

    """

    def __init__(self, pdb_file):
        """Constructor

        Arguments:
            pdb_file -- the location of the PDB file of the protein.
        """

        # Check if the pdb file is gzipped, gunzip it to a temp file if it is.
        if os.path.splitext(pdb_file)[1] == ".gz":
            (tmp_f, tmp_fn) = tempfile.mkstemp(".pdb")
            os.system("gunzip -c " + pdb_file + " > " + tmp_fn)
            pdb_file = tmp_fn

        self.__structure = Structure(pdb_file)
        self.pdb_id = self.__structure.pdb_code

        # Need to check if chains are incorrectly split like for 1ptm.
        chain_locs = {}
        tmp_chains = []
        chain_ids = []
        prev_chain_id = None

        _chains = self.__structure.peptide_chains

        # BUG 2014-12-30:
        # MK reported failure on sample CASP PDB file (SampleCASPFile.pdf),
        # which doesn't contain any chain identifiers.
        #
        # SOL:
        # Scan for chains w/ missing chain identifiers, assign placeholder chain
        # identifier '?' (a single question mark) to these chains, of course,
        # this assumes that such an error will occur only for PDB files with a
        # single chain. Not yet sure how to appropriately handle multi-chain PDB
        # files with multiple missing chain identifiers.
        
        for chain in self.__structure.peptide_chains:
            chain.chain_id = default_chain_id(chain.chain_id)

        # Fix broken up chains.
        if len(_chains) > 1:
            for i in xrange(len(_chains)):
                while _chains[i+1].chain_id == _chains[i].chain_id:
                    self.__structure.joinPeptideChains(i, i+1)

        # Create a wrapper object around the ScientificPython Chain objects.
        self.chains = [Chain(chain) for chain in \
                       self.__structure.peptide_chains]

        self.__parse_dssp_stride(pdb_file)
        self.__calc_strand_directions()

        # Delete the temp file if it was used.
        try:
            os.close(tmp_f)
            os.system("\\rm -f " + tmp_fn)
        except NameError:
            pass

        # If no residues were parsed by the ScientificPython parser then exit
        # failure. This can occur because pdb files may only contain nucleic
        # acids.
        if len(self.get_polypeptide()) == 0:
            sys.stderr.write("error: " + self.pdb_id + \
                             " contains no amino acids\n")
            sys.exit(1)

    def __calc_strand_directions(self):
        """Calculate the direction of the beta-strands. The direction of a
        beta-strand is taken as the vector implied by the first residue to the
        last residue of the strand. This direction vector is then assigned to
        each residue of the strand.

        """

        prev_sse = None
        pp = self.get_polypeptide()

        for i in xrange(len(pp)):
            if pp[i].stride_sse == "E" and prev_sse != "E":
                # Reached the start of a new strand.
                first_v = pp[i].get_atom("CA").position
                first_i = i
            elif pp[i].stride_sse != "E" and prev_sse == "E":
                # Reached the end of a strand.
                last_v = pp[i-1].get_atom("CA").position
                strand_direction = last_v - first_v

                # If the strand length is only one residue then use the
                # direction vector implied by the N to C teriminus.
                if first_i == i - 1:
                    strand_direction = pp[first_i].get_atom("C").position - \
                            pp[first_i].get_atom("N").position

                # Assign the direction to each residue in the strand.
                for j in xrange(first_i, i):
                    pp[j].strand_direction = strand_direction

            prev_sse = pp[i].stride_sse

        return

    def __parse_dssp_stride(self, pdb_file):
        """Parse output of dssp and stride.

        Arguments:
            pdb_file -- the location of the PDB file of the protein.
        """

        # Put contents of dssp output into a temp file.
        (tmp_f, tmp_fn) = tempfile.mkstemp(".stridedssp")
        #os.system("dsspcmbi -w " + pdb_file + " " + tmp_fn) 
        os.system("dssp " + pdb_file + " " + tmp_fn) 
        output = open(tmp_fn, "r")
        self.__dssp_lookup = {}

        # Skip to the important records.
        for line in output:
            if line[2] != "#":
                continue
            break

        # Parse the contents of the temp file.
        for line in output:

            # Skip over chain terminating characters.
            if line[13] == "!":
                continue

            dssp_id = line[:5].strip()
            pdb_res_id = line[5:10].strip()
            chain_id = default_chain_id(line[11])
            dssp_sse = line[16]
            bp1 = line[25:31].strip()
            bp2 = line[31:37].strip()

            # Check if the Scientific.IO.PDB.Structure parser recognised the
            # residue.
            try:
                res = self.get_pdb_residue(chain_id, pdb_res_id)
            except (AttributeError, KeyError):
                # Residue not recognised, move on.
                print traceback.format_exc()
                continue

            # Check if dssp doesn't assign an sse.
            if dssp_sse.isspace():
                dssp_sse = None

            # Check for empty bridge partners.
            if bp1 == "0":
                bp1 = None
            if bp2 == "0":
                bp2 = None

            # Assign the dssp-parsed values to the residue, residue is
            # modified in-place.
            res.dssp_id = dssp_id
            res.dssp_sse = dssp_sse

            res.bp1 = bp1
            res.bp2 = bp2

            self.__dssp_lookup[dssp_id] = res

        output.close()

        # Parse stride output.
        # Hydrogen bond partners are determined here.
        os.system("stride -h " + pdb_file + " > " + tmp_fn)
        output = open(tmp_fn, "r")

        for line in output:
            rec_type = line[:3]
            
            if rec_type == "ACC" or rec_type == "DNR":
                hbp1 = line[11:15].strip()
                hbp2 = line[31:35].strip()

                hbp1_chain_id = default_chain_id(line[9])
                hbp2_chain_id = default_chain_id(line[29])

                assert(hbp1_chain_id.strip())
                assert(hbp2_chain_id.strip())

                #print hbp1_chain_id, hbp2_chain_id
                 
                # Check if the current residue has been parsed.
                try:
                    curr_res = self.get_pdb_residue(hbp1_chain_id, hbp1)
                except (AttributeError, KeyError):
                    # Residue not parsed, move on.
                    continue

                # Check if the partner residue has been parsed.
                try:
                    hbond_res = self.get_pdb_residue(hbp2_chain_id, hbp2)
                except (AttributeError, KeyError):
                    # Residue not parsed, move on.
                    continue

                # Check if the partner residue was parsed by dssp.
                if hbond_res.dssp_id:
                    curr_res.hbonds.append(hbond_res.dssp_id)
                if curr_res.dssp_id:
                    hbond_res.hbonds.append(curr_res.dssp_id)

            if rec_type == "ASG":
                stride_sse = line[24]
                chain_id = default_chain_id(line[9])
                pdb_res_id = line[11:15].strip()

                # Check if the current residue was parsed by dssp.
                try:
                    curr_res = self.get_pdb_residue(chain_id, pdb_res_id) 
                    curr_res.stride_sse = stride_sse
                except (AttributeError, KeyError):
                    # Residue not parsed, move on.
                    pass

        output.close()

        # Delete temp file.
        os.close(tmp_f)
        os.system("\\rm -f " + tmp_fn)

        return

    def __iter__(self):
        """Return an iterator of the chains in a protein.

        """

        return iter(self.chains)

    def get_chain(self, chain_id):
        """Get a chain by its chain id.

        Arguments:
            chain_id -- the chain's ID (e.g. "A")

        """

        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain

    def get_dssp_residue(self, dssp_id):
        """Get a Residue object by its dssp_id.

        Arguments:
            dssp_id -- the ID number assigned by DSSP.

        """

        return self.__dssp_lookup[dssp_id] 

    def get_pdb_residue(self, chain_id, pdb_res_id):
        """Get a Residue object by its pdb residue id, because pdb residue
        ids aren't globally unique, you must specify the chain id.

        Arguments:
            chain_id -- the ID of the chain that the residue is located in.
            pdb_res_id -- the PDB residue ID of the residue from the PDB file.

        """

        return self.get_chain(chain_id).get_pdb_residue(pdb_res_id)

    def get_polypeptide(self, chains=[]):
        """Return a flat list of all the residues in the protein. Sometimes it
        is easier to work with a flat list.

        Arguments:
            chains -- a list of chain IDs of residues to be included in the
            polypeptide.

        """
    
        # If no chain was specified, return all the residues from all chains.
        if chains == []:
            return [residue for chain in self.chains for residue in \
                    chain.residues]

        polypeptide = []

        for chain in self.chains:
            if chain.chain_id not in chains:
                continue
            for residue in chain.residues:
                polypeptide += [residue]

        return polypeptide


class Chain(object):
    """A protein chain.

    """

    def __init__(self, peptide_chain):
        """Constructor.

        Arguments:
            peptide_chain -- the peptide chain object obtained when calling
            Scientific.IO.PDB.Structure.peptide_chains.

        """

        self.__peptide_chain = peptide_chain
        self.residues = [Residue(residue, peptide_chain.chain_id) \
                         for residue in peptide_chain]
        self.chain_id = peptide_chain.chain_id
        self.sequence = peptide_chain.sequence

        self.__init_pdb_lookup()
        return

    def __iter__(self):
        return iter(self.residues)

    def __str__(self):
        return str(self.__peptide_chain)

    def __init_pdb_lookup(self):
        self.__pdb_lookup = {}

        for residue in self.residues:
            self.__pdb_lookup[str(residue.number)] = residue

        return

    def get_pdb_residue(self, pdb_res_id):
        """Get a residue by its pdb residue number.

        Arguments:
            pdb_res_id -- the residue ID of the residue in the PDB file.
        
        """
        return self.__pdb_lookup[str(pdb_res_id)]


class Residue(object):
    """An amino acid residue.
    """

    def __init__(self, residue, chain_id):
        """Constructor.

        Arguments:
            residue -- the Scientific.IO.PDB.AminoAcidResidue object obtained
            from ScientificPython's handling of PDB parsing.
            chain_id -- the chain ID of the residue.

        """

        # Inherited attributes.
        self.__residue = residue
        self.name = residue.name
        self.number = residue.number
        self.atoms = residue.atoms
        self.atom_list = residue.atom_list
        self.chain_id = chain_id

        # New attributes.
        self.dssp_id = None
        self.dssp_sse = None
        self.stride_sse = None
        self.strand_direction = None
        self.bp1 = None
        self.bp2 = None
        self.hbonds = []

        return

    def __str__(self):
        return str(self.__residue)

    def get_atom(self, atom_name):
        """Return the specified atom, or the next closest one if it is
        missing.

        Arguments:
            atom_name -- the name of the atom in the residue.
        
        """

        try:
            return self.atoms[atom_name]
        except KeyError:
            pass

        # If the specified atom isn't there, then use the next best one.
        for atom_name in representative_atoms:
            try:
                return self.atoms[atom_name]
            except KeyError:
                continue


class Matrix(object):
    """Wrapper class around numpy.matrix that allows cartesian-coordinate
    indexing.

    """

    def __init__(self, data):
        """Constructor.

        Arguments:
            data -- a two-dimensional array or list to initialise the matrix.

        """

        self.__matrix = matrix(data)

    def __len__(self):
        return len(self.__matrix)
    
    def get(self, x, y):
        """Get a matrix value using a cartesian-coordinate system.

        Arguments:
            x -- the x coordinate of the matrix.
            y -- the y coordinate of the matrix.

        Returns:
            the value located at the (x, y) coordinate of the matrix.

        """

        if x > self.__matrix.shape[1] - 1 or y > self.__matrix.shape[0] - 1:
	        raise IndexError

        return self.__matrix[len(self.__matrix) - y - 1, x]
    
    def set(self, x, y, val):
        """Set a matrix value at a particular coordinate.

        Arguments:
            x -- the x coordinate of the matrix.
            y -- the y coordinate of the matrix.
            val -- the value to be set.

        """

        if x > self.__matrix.shape[1] - 1 or y > self.__matrix.shape[0] - 1:
            raise IndexError

        self.__matrix[len(self.__matrix) - y - 1, x] = val

    def get_matrix(self):
        """Get the numpy.matrix object that this class wraps around.

        Returns:
            the numpy.matrix object that this class wraps.

        """

        return self.__matrix

    def save_matrix(self, fn):
        """Save the contents of the contact matrix in ASCII format to a file.

        Arguments:
            fn -- the filename to save the matrix to.

        """

        f = open(fn, "w")

        for i in xrange(len(self)):
            line = ""
            for j in xrange(len(self)):
                line += str(self.get(i, j)) + " "
            f.write(line + "\n")
        f.close()

        return


class DistanceMatrix(Matrix):
    """Subclass of Matrix is specialised for representing real-valued distance
    maps. A DistanceMatrix is a N x N matrix.

    """

    def __init__(self, pp, metric, seq_separation):
        """Constructor.

        Arguments:
            pp -- the list of residues of the protein.
            metric -- the inter-residue distance metric to be used.
            seq_separation -- the minimum sequence separation of which a
                contact is defined.

        """

        # Call superclass constructor with a blank matrix of floats.
        Matrix.__init__(self, zeros((len(pp), len(pp)), float))

        self.pp = pp
        self.metric = metric
        self.seq_separation = seq_separation

        # Calculate the distance matrix.
        self.__calc_distance_matrix()

        return

    def __calc_distance_matrix(self):
        """Calculate the distance matrix.

        """

        # Populate the matrix.
        for x in xrange(len(self.pp)):
            for y in xrange(len(self.pp)):
                if abs(x - y) >= self.seq_separation:
                    self.set(x, y, calc_distance(self.pp[x], self.pp[y],
                                                 self.metric))

        return

    def plot(self, plot_title=None, colourbar=False, chain_boundaries=False,
             noticks=False, output=None):
        """Plot the distance map.

        Arguments:
            plot_title -- the title that will be displayed at the top of the
                map.
            colourbar -- if set to True, then a colour bar will be drawn on
                the right-hand side of the map.
            chain_boundaries -- if set to True, then vertical lines indicating
                the locations of chains will be drawn on the map.
            noticks -- if set to True, then the sequence index numbers will
                not be drawn on the axes.
            output -- the file for the plot to be saved to.

        """

        dmatrix = self.get_matrix().tolist()
        dmatrix.reverse()

        draw_axes()

        im = pylab.imshow(dmatrix, interpolation="nearest", origin="lower")
        im.figure.set_size_inches(12, 12, forward=True)

        draw_grid()
        draw_title(plot_title)
        draw_chain_boundaries(self.pp, chain_boundaries)
        draw_ticks(noticks)

        if colourbar:
            cb = pylab.colorbar(shrink=0.8)

            # Change the font family of colour bar tick labels.
            for tick_labels in cb.ax.get_yticklabels():
                tick_labels.set_name(FONT_NAME)

        if output:
            pylab.savefig(output)
        else:
            pylab.show()

        return


class ContactMatrix(Matrix):
    """Subclass of Matrix that is specialised for representing binary-valued
    contact maps. A ContactMatrix is a N x N square matrix.

    """

    def __init__(self, pp, metric, threshold, min_threshold, seq_separation):
        """Constructor.

        Arguments:
            pp -- the list of residues of the protein.
            metric -- the inter-residue distance metric to be used.
            min_threshold -- the minimum distance threshold (in Angstroms) of
                which a contact is defined.
            seq_separation -- the minimum sequence separation of which a
                contact is defined.

        """

        # Call superclass constructor with a blank matrix of ints.
        Matrix.__init__(self, zeros((len(pp), len(pp)), int)) 

        # Attributes.
        self.pp = pp
        self.metric = metric
        self.threshold = threshold
        self.min_threshold = min_threshold
        self.seq_separation = seq_separation

        # Calculate the contact matrix.
        self.__calc_contact_matrix()

        return

    def __calc_contact_matrix(self):
        """Calculate the contact matrix.

        """

        # Populate the matrix.
        x = 0
        n = len(self.pp)

        for x in xrange(len(self.pp)):
            for y in xrange(len(self.pp)):
                dist = calc_distance(self.pp[x], self.pp[y], self.metric)

                if dist < self.threshold and dist >= self.min_threshold and \
                abs(x - y) >= self.seq_separation:
                    self.set(x, y, 1)

        return


    def print_contact_list(self, fh):
        """ Print the contact list
        
        Write the contact list to fh in the following format:

        n

        ri<tab>rj
        ...

        where n is the number of residues and each line has ri
        and rj, residue numbers in contact separate by tab char.

        Arguments:
            fh - open for write filehandle to write to
        """
        fh.write(str(len(self)) + '\n\n')

        for i in xrange(len(self)):
            for j in xrange(i+1, len(self)):
                if self.get(i,j) > 0:
                    fh.write(str(i) + '\t' + str(j) + '\n')

    def plot(self, sse=None, hbonds=False, regions=None, chain_boundaries=None,
             plot_title=None, noticks=False, nolegend=False, output=None):
        """Plot the contact matrix.

        Arguments:
            sse -- if set to "DSSP" or "STRIDE", the contact map will be
                annotated with secondary structure information using the
                respective algorithm.
            hbonds -- if set to True, hydrogen-bond information between
                contacts will be drawn on the map.
            regions -- a string specifying which regions of residues will be
                highlighted on the map.
            plot_title -- the title that will be displayed at the top of the
                map.
            chain_boundaries -- if set to True, then vertical lines indicating
                the locations of chains will be drawn on the map.
            plot_title -- the title that will be displayed at the top of the
                map.
            noticks -- if set to True, then the sequence index numbers will
                not be drawn on the axes.
            nolegend -- if set to True, sse colour code legend will not be
                displayed on the map.
            output -- the file for the plot to be saved to.

        """
        
        # Regular contact coordinates.
        xs = []
        ys = []

        # Hydrogen bond coordinates.
        hbonds_xs = []
        hbonds_ys = []

        draw_figure()
        draw_axes()
        draw_grid()

        # Try to scale the markersize according to the length of the protein.
        ms = 700.0 / len(self.pp)

        sse_locs = {"H" : [], "B" : [], "E" : [], "G" : [], "I" : [],
                    "T" : [], "S" : [], "C" : []}

        sse_labels = {"H" : "A-helix", "B" : "Isolated B-bridge", \
                      "E" : "B-strand", "G" : "3/10-helix", \
                      "I" : "Pi-helix", "T" : "Turn", \
                      "S" : "Bend", "C" : "Coil"}

        for i, res_one in enumerate(self.pp):
            for j, res_two in enumerate(self.pp):
                if (i == j and not sse) or abs(i-j) < self.seq_separation:
                    continue

                # Obtain the sse annotations.
                if sse == "DSSP":
                    sse_i = self.pp[i].dssp_sse
                    sse_j = self.pp[j].dssp_sse
                else:
                    sse_i = self.pp[i].stride_sse
                    sse_j = self.pp[j].stride_sse

                if self.get(i, j) > 0:
                    xs += [i+1]
                    ys += [j+1]

                    if sse and sse_i and sse_i == sse_j:
                        sse_locs[sse_i] += [(i+1, j+1)]

                    if hbonds and self.pp[j].dssp_id in self.pp[i].hbonds:
                        hbonds_xs += [i+1]
                        hbonds_ys += [j+1]

                elif i == j and sse:
                    if sse_i and sse_i == sse_j:
                        sse_locs[sse_i] += [(i+1, j+1)]

        pylab.plot(xs, ys, "ks", ms=ms,  mew=0, mfc="black", zorder=1, aa=False)

        # Plot the secondary structure annotations if chosen.
        if sse:
            for key in sse_locs.keys():
                locs = sse_locs[key.upper()]
                sse_xs = [x for (x, _) in locs]
                sse_ys = [y for (_, y) in locs]

                if len(locs) == 0:
                    continue

                pylab.plot(sse_xs, sse_ys, mfc=SSE_COLOUR[key], marker="s",
                           mew=0, mec=SSE_COLOUR[key], ms=ms, aa=False,
                           label=sse_labels[key], lw=0, zorder=5)

        # Plot hbonds if specified.
        if hbonds:
            pylab.plot(hbonds_xs, hbonds_ys, "ro", ms=3, aa=False, zorder=11,
                      label="Hydrogen-bond")

        draw_title(plot_title)
        draw_highlights(regions)
        draw_chain_boundaries(self.pp, chain_boundaries)
        draw_ticks(noticks)

        # Draw a legend if specified.
        if sse and not nolegend:
            pylab.legend(loc=0, markerscale=1.0, handletextsep=0, numpoints=1)

        if output:
            pylab.savefig(output)
        else:
            pylab.show()

        return


class HydrogenBondMatrix(Matrix):
    """Subclass of Matrix that is specialised for representing hydrogen-bonds
    between residues.

    """

    def __init__(self, pp):
        """Constructor.

        Arguments:
            pp -- the list of residues of the protein.

        """

        # Call superclass constructor with a blank matrix of ints.
        Matrix.__init__(self, zeros((len(pp), len(pp)), int))

        self.pp = pp
        self.__calc_hydrogen_bond_matrix()

        return

    def __calc_hydrogen_bond_matrix(self):
        """Calculate the hydrogen bond matrix.

        """

        for i in xrange(len(self.pp)):
            for j in xrange(len(self.pp)):
                if self.pp[j].dssp_id in self.pp[i].hbonds:
                    self.set(i, j, 1)

        return

    def plot(self, plot_title=None, chain_boundaries=None,
             regions=None, noticks=False, output=None):
        """Plot the hydrogen bond plot.

        Arguments:
            regions -- a string specifying which regions of residues will be
                highlighted on the map.
            plot_title -- the title that will be displayed at the top of the
                map.
            chain_boundaries -- if set to True, then vertical lines indicating
                the locations of chains will be drawn on the map.
            plot_title -- the title that will be displayed at the top of the
                map.
            noticks -- if set to True, then the sequence index numbers will
                not be drawn on the axes.
            output -- the file for the plot to be saved to.

        """

        xs = []
        ys = []

        for i in xrange(len(self.pp)):
            for j in xrange(len(self.pp)):
                if i == j:
                    continue

                if self.get(i, j) > 0:
                    xs += [i+1]
                    ys += [j+1]

        draw_figure()
        draw_axes()
        draw_grid()

        ms = 700.0 / len(self.pp)

        pylab.plot(xs, ys, "rs", ms=ms, mew=0.0, aa=False, zorder=10)

        draw_highlights(regions)
        draw_title(plot_title)
        draw_chain_boundaries(self.pp, chain_boundaries)
        draw_ticks(noticks)

        if output:
            pylab.savefig(output)
        else:
            pylab.show()

        return


def calc_angle(v1, v2):
    """Calculate the angle (in degrees) between two three-dimensional
    direction vectors.

    Arguments:
        v1 -- the first vector.
        v2 -- the second vector.

    Returns:
        the angle between v1 and v2 in degrees.

    """

    v1 = array(v1)
    v2 = array(v2)
    v1 /= norm(v1)
    v2 /= norm(v2)

    if v1.tolist() == v2.tolist():
        return 0.0

    return arccos(dot(v1, v2) / (norm(v1)*norm(v2))) * 180 / pi


def calc_eucl_distance(v1, v2):
    """Calculate the Euclidean distance between two vectors.

    Arguments:
        v1 -- the first vector.
        v2 -- the second vector.

    Returns:
        the Euclidean distance between v1 and v2.

    """
    
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + 
                      (v1[2] - v2[2]) ** 2)


def calc_COM(res):
    """Calculate the center of mass of the entire residue.

    Arguments:
        res -- the residue of interest.

    Returns:
        the coordinates of the all-atom center of mass.

    """

    com_coord = array([0.0, 0.0, 0.0])
    res_mass = 0.0

    for atom in res.atoms.keys():
        com_coord += MASS[atom[0]] * res.atoms[atom].position
        res_mass += MASS[atom[0]]

    if res_mass == 0:
        return res.get_atom("CA").position

    return com_coord / res_mass


def calc_RCOM(res):
    """Calculate the center of mass of only the side chain (R-group) residues,
    this excludes the CA atom.

    Arguments:
        res -- the residue of interest.

    Returns:
        the coordinates of the side chain center of mass.

    """

    # Glycine doesn't have a side chain, use the CA atom instead.
    if res.name == "GLY":
        return res.get_atom("CA").position

    side_coord = array([0.0, 0.0, 0.0])
    side_mass = 0.0

    # Select only the side chain atoms.
    ignored_atoms = ["N", "CA", "C", "O", "H"]
    noticed_atoms = [atom for atom in res.atoms.keys() if atom not in
                     ignored_atoms]

    for atom in noticed_atoms:
        side_coord += MASS[atom[0]] * res.atoms[atom].position
        side_mass += MASS[atom[0]]

    if side_mass == 0:
        return res.get_atom("CA").position

    return side_coord / side_mass


def calc_minimum_distance(res_x, res_y):
    """Calculate the minimum distance between two residues using an all-atom
    representation.

    Arguments:
        res_x -- the first residue.
        res_y -- the second residue.

    Returns:
        the minimal distance between the two residues.

    """

    distances = []

    for atom_x in res_x.atoms.values():
        for atom_y in res_y.atoms.values():
            distances += [calc_eucl_distance(atom_x.position, atom_y.position)]

    distances.sort()

    return distances[0]


def calc_min_VDW_distance(res_x, res_y):
    """Calculate the minimum VDW distance between two residues.

    Arguments:
        res_x -- the first residue.
        res_y -- the second residue.

    Returns:
        the minimum VDW distance between the two residues.

    """

    distances = []

    for atom_x in res_x.atoms.values():
        for atom_y in res_y.atoms.values():
            distances += [calc_eucl_distance(atom_x.position, atom_y.position) -
                          VDW_RADII[atom_x.name[0]] -
                          VDW_RADII[atom_y.name[0]]]

    distances.sort()

    return distances[0]


def calc_distance(res_x, res_y, metric):
    """Calculate the distance between two residues according to specified
    metric.

    Arguments:
        res_x -- the first residue.
        res_y -- the second residue.
        metric -- the inter-residue distance metric to use:
            "CA" --> ca-ca distance.
            "CB" --> cb-cb distance.
            "COM" --> distances between residue centres of mass.
            "RCOM" --> distances between residue side-chain centres of mass.
            "MIN" --> minimum distance between residue atoms.
            "VDW" --> minimum distance between residue atoms' VDW radii.

    Returns:
        the distance between two residues according to a specified
        inter-residue distance metric.

    """

    if metric == "CA":
        return calc_eucl_distance(res_x.get_atom("CA").position,
                                  res_y.get_atom("CA").position)
    elif metric == "CB":
        return calc_eucl_distance(res_x.get_atom("CB").position,
                                  res_y.get_atom("CB").position)
    elif metric == "COM":
        return calc_eucl_distance(calc_COM(res_x), calc_COM(res_y))

    elif metric == "RCOM":
        return calc_eucl_distance(calc_RCOM(res_x), calc_RCOM(res_y))

    elif metric == "MIN":
        return calc_minimum_distance(res_x, res_y)

    elif metric == "VDW":
        return calc_min_VDW_distance(res_x, res_y)


def draw_axes():
    """Draw the axis on the plot with appropriately sized margins.

    """

    pylab.axes([0.05, 0.05, 0.9, 0.9])

    return


def draw_figure():
    """Draw an appropriately sized empty figure.

    """

    pylab.figure(figsize=(12,12))

    return


def draw_grid():
    """Draw a grid on the plot.

    """

    pylab.grid(alpha=0.5)

    return


def draw_highlights(regions):
    """Highlight specified residue region(s) on the plot.

    Arguments:

        regions -- the specific residues that will be highlighted on the plot
        according to the following format examples:

            "33:48" for residues 33 to 48
            "33:" for residues 33 onwards
            "33:48,67:78" for residues 33 to 48 and 67 to 78

    """

    if not regions:
        return

    regions = regions.split(",")
    ax = pylab.gca()

    ylim = ax.get_ylim()
    xlim = ax.get_xlim()

    for r in regions:
        start = r.split(":")[0]
        end = r.split(":")[1]

        if start == '':
            start = xlim[0]

        if end == '':
            end = xlim[1]

        start = int(start)
        end = int(end)
        
        if start < xlim[0]:
            start = xlim[0]
        elif start > xlim[1]:
            start = xlim[1]

        if end > xlim[1]:
            end = xlim[1]
        elif end < xlim[0]:
            end = xlim[0]

        xs = [start, end, end, start]
        ys = [ylim[0], ylim[0], ylim[1], ylim[1]]

        pylab.fill(xs, ys, lw=0, fc="#eaeaea")
        pylab.fill(ys, xs, lw=0, fc="#eaeaea")

    return


def draw_chain_boundaries(pp, chain_boundaries):
    """Plot the chain boundaries of multi-chain proteins.

    Arguments:
        pp -- flat list of residues in the plot.
        chain_boundaries -- if set to True then draw the chain boundaries.

    """

    if not chain_boundaries:
        return

    prev_chain = ""

    for i, res in enumerate(pp):
        if res.chain_id != prev_chain:
            pylab.axvline(x=i, c="black", lw=0.5, zorder=100)
            pylab.axhline(y=i, c="black", lw=0.5, zorder=100)
        prev_chain = res.chain_id

    return


def draw_title(plot_title):
    """Draw the plot title at the top of the plot.

    Arguments:
        plot_title -- the plot title.

    """

    if plot_title:
        pylab.title(plot_title, name=FONT_NAME, fontweight="bold",
                    fontsize="xx-large", verticalalignment="bottom")

    return


def draw_ticks(noticks=False):
    """Draw the residue indices on the axes.

    Arguments:
        noticks -- if set to False then draw the residue indices on the axes.

    """

    if not noticks:
        pylab.xticks(name=FONT_NAME, fontsize="xx-large")
        pylab.yticks(name=FONT_NAME, fontsize="xx-large")
    else:
        pylab.setp(pylab.gca(), xticklabels=[], yticklabels=[])

    return


def make_plot():
    """Use plot.py as a standalone script to plot contact maps from the
    command line.

    """

    # Specify the options.
    parser = OptionParser("%prog [ --cmap | --dmap | --hbplot ] --pdb <PDB file> [options]" )

    parser.add_option("--hbplot", dest="hbplot", default=False,
                      action="store_true",
                      help="create a hydrogen-bond plot")

    parser.add_option("--dmap", dest="dmap", default=False,
                      action="store_true",
                      help="create a distance map")

    parser.add_option("--cmap", dest="cmap", default=True,
                      action="store_true",
                      help="create a contact map")

    parser.add_option("--pdb", dest="pdb",
                      help="pdb file of the protein (e.g. 1UBQ.pdb) " +
                      "(required).")

    parser.add_option("--output", dest="output", default=None,
                      help="file to save the plot to (e.g. *.eps, " +
                      "*.pdb, *.png, *.ps, or *.svg), if not specified, " +
                      "the plot will be displayed on the screen in " + 
                      "interactive mode.")

    parser.add_option("--chains", dest="chains", default="all",
                      help="selected chains to plot (e.g. --chains ABC), " +
                      "(default=all)")

    parser.add_option("--threshold", dest="threshold", default=8.0,
                      type="float", help="contact threshold (default=8.0).")

    parser.add_option("--metric", dest="metric", default="CA",
                      choices=["CA", "CB", "COM", "RCOM", "MIN", "VDW"],
                      help="distance metric. Choices: " +
                      "CA - alpha carbons, " +
                      "CB - beta carbons (CASP standard), " +
                      "COM - residue centers of mass, " +
                      "RCOM - side chain centers of mass, " +
                      "MIN - minimum atom distance, " +
                      "VDW - minimum van der Waals distance, (default=CA).")

    parser.add_option("--min_threshold", dest="min_threshold", default=0.0,
                      type="float",
                      help="minimum distance threshold (default=0.0).")

    parser.add_option("--seq_separation", dest="seq_separation", default=0,
                      type="int",
                      help="minimum sequence separation between residues " +
                      "(default=0).")

    parser.add_option("--highlight", dest="regions", default=None,
                      help="highlight specific residues " + 
                      "(disabled by default)")

    parser.add_option("--sse", dest="sse", default=None,
                      choices=["STRIDE", "DSSP"],
                      help="annotate contacts with SSE information " +
                      "(disabled by default). Choices: STRIDE or DSSP.")

    parser.add_option("--hbonds", dest="hbonds", action="store_true",
                      default=False,
                      help="plots hydrogen bonds if enabled " + 
                      "(disabled by default).")

    parser.add_option("--title", dest="plot_title", default=None,
                      help="the plot title.")

    parser.add_option("--colourbar", dest="colourbar", default=False,
                      action="store_true",
                      help="display colour bar for distance maps " +
                      "(disabled by default).")

    parser.add_option("--chain_boundaries", dest="chain_boundaries",
                      default=None, action="store_true",
                      help="draw chain boundaries. (disabled by default)")

    parser.add_option("--noticks", dest="noticks", default=False,
                      action="store_true",
                      help="turn off tick labelling on the axes " +
                      "(ticks are displayed by default).")

    parser.add_option("--ascii", dest="ascii", default=False,
                      action="store_true",
                      help="save the plot as an ascii file.")

    parser.add_option("--cmaplist", dest="cmaplist", default=False,
                      action="store_true",
                      help=" save the contact map as an ASCII list of" 
                           " contact residue pairs.")
                     

    parser.add_option("--nolegend", dest="nolegend", default=False,
                      action="store_true",
                      help="turn off the sse colour code legend " +
                      "(legend is displayed by default if --sse option is" +
                      "specified")
    
    # Parse the options.
    (options, _) = parser.parse_args()

    if not options.pdb:
        sys.stderr.write('A PDB file must be specified with --pdb\n')
        parser.print_help()
        sys.exit(1)

    # Parse chains.
    if options.chains == "all":
        chains = []
    else:
        chains = [c.upper() for c in options.chains]

    protein = Protein(options.pdb)
    pp = protein.get_polypeptide(chains)

    if options.hbplot:
        hbmatrix = HydrogenBondMatrix(pp)

        if options.ascii:
            hbmatrix.save_matrix(options.output)
        else:
            hbmatrix.plot(regions=options.regions, plot_title=options.plot_title,
                          chain_boundaries=options.chain_boundaries,
                          noticks=options.noticks, output=options.output)
    elif options.dmap:
        dmatrix = DistanceMatrix(pp, metric=options.metric,
                                 seq_separation=options.seq_separation)

        if options.ascii:
            dmatrix.save_matrix(options.output)
        else:
            dmatrix.plot(plot_title=options.plot_title,
                         colourbar=options.colourbar,
                         chain_boundaries=options.chain_boundaries,
                         noticks=options.noticks,
                         output=options.output)

    elif options.cmap:

        if not options.output:
            sys.stderr.write("ERROR: please specify an output file using"
                             " --ouput <filename>.\n")
            sys.exit(-1)

        cmatrix = ContactMatrix(pp, metric=options.metric,
                                threshold=options.threshold,
                                min_threshold=options.min_threshold,
                                seq_separation=options.seq_separation)

        if options.ascii:
            cmatrix.save_matrix(options.output)
        elif options.cmaplist:
            if options.output:
                fh = open(options.output, 'w')
            else:
                fh = sys.stdout
            cmatrix.print_contact_list(fh)
            fh.close()
        else:
            cmatrix.plot(sse=options.sse, hbonds=options.hbonds,
                         regions=options.regions,
                         chain_boundaries=options.chain_boundaries,
                         plot_title=options.plot_title,
                         noticks=options.noticks,
                         nolegend=options.nolegend,
                         output=options.output)

if __name__ == "__main__":
    make_plot()

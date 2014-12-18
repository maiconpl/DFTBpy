"""
This file is part of DFTBpy.

    DFTBpy is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DFTBpy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DFTBpy. If not, see <http://www.gnu.org/licenses/>.
"""

from math import sqrt
from model.slater_koster_transformations import SlaterKosterTransformations
import sys

class DFTBMatrix(object):

    def __init__(self, dist_x,
                 dist_y,
                 dist_z,
                 dist_matrix,
                 atom_symb,
                 atom_numb,
                 atom_type_numb,
                 max_all_momentum,
                 atom_orb_sym,
                 total_orb,
                 H_diag_dftb,
                 dist_list,
                 dist_step_list,
                 H_list, 
                 S_list,
                 H_cutoff_list,
                 debug):
        
        self.dist_x = dist_x
        self.dist_y = dist_y
        self.dist_z = dist_z
        self.dist_matrix = dist_matrix
        self.atom_numb = atom_numb
        self.atom_symb = atom_symb
        self.atom_type_numb = atom_type_numb
        self.max_all_momentum = max_all_momentum
        self.atom_orb_sym = atom_orb_sym
        self.total_orb = total_orb
        self.H_diag_dftb = H_diag_dftb
        self.H_cutoff_list = H_cutoff_list
        self.debug = debug
        
        self.H_dftb = []
        self.S_dftb = []
        
        self.sk_transformations = SlaterKosterTransformations(dist_list,
                                                              dist_step_list,
                                                              H_list,
                                                              S_list,
                                                              max_all_momentum)

    def build_H_S_matrix(self):

        self.H_dftb = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)] # defining the dimension of the matrix.
        self.S_dftb = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)] # defining the dimension of the matrix.

        # These loops contemplates only the creation of the right (upper) triangle H and S matrix.
        for itotal_orb in xrange(self.total_orb):
            for jtotal_orb in xrange(itotal_orb, self.total_orb):

                dist_x = self.dist_x[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]]
                dist_y = self.dist_y[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]]
                dist_z = self.dist_z[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]]

                # Condiction of same orbital in the same center (diagonal of the matrix):
                if itotal_orb == jtotal_orb:
                   self.H_dftb[itotal_orb][itotal_orb] = self.H_diag_dftb[itotal_orb]
                   self.S_dftb[itotal_orb][itotal_orb] = 1.0

                # Condition of different orbitals in the same center:
                if itotal_orb != jtotal_orb and self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]] == 0.0:
                   self.H_dftb[itotal_orb][jtotal_orb] = 0.0
                   self.S_dftb[itotal_orb][jtotal_orb] = 0.0

                # Condition when the interatomic distance(s) is/are below the integral cutoff for <s|H|s>: usually below 0.5 bohr.
                if itotal_orb != jtotal_orb and self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]] != 0.0 and self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]] < self.H_cutoff_list[self.atom_type_numb[itotal_orb]][self.atom_type_numb[jtotal_orb]][9][0]:

                   print ""
                   print "Program stopped."
                   print "Warning! The distance between an/some atom(s) in the system is it too small, below " + str(self.H_cutoff_list[self.atom_type_numb[itotal_orb]][self.atom_type_numb[jtotal_orb]][9][0])  + " bohr."
                   print "The distance between the atoms", str(self.atom_numb[itotal_orb] + 1) + " and " + str(self.atom_numb[jtotal_orb] + 1) + " is:", self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]], "bohr."
                   print "Please! Check your coordinates."
                   print ""
                   sys.exit()

                # Condition when the interatomic distance(s) is/are greater than that/those supported by the SLAKO: <s|H|s>.
                if itotal_orb != jtotal_orb and self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]] > self.H_cutoff_list[self.atom_type_numb[itotal_orb]][self.atom_type_numb[jtotal_orb]][9][1]:

                   self.H_dftb[itotal_orb][jtotal_orb] = 0.0
                   self.S_dftb[itotal_orb][jtotal_orb] = 0.0

                # Condition of the interaction of atomic orbitals of one center with another center:
                if itotal_orb != jtotal_orb and self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]] != 0.0:

                   dist_atoms = self.dist_matrix[self.atom_numb[itotal_orb]][self.atom_numb[jtotal_orb]]
                   integral_symmetry = self.atom_orb_sym[itotal_orb] + self.atom_orb_sym[jtotal_orb]
                   self.H_dftb[itotal_orb][jtotal_orb], self.S_dftb[itotal_orb][jtotal_orb] = self.sk_transformations.get_skf(self.atom_type_numb[itotal_orb],
                                                                                                                               self.atom_type_numb[jtotal_orb],
                                                                                                                               integral_symmetry,
                                                                                                                               dist_atoms,
                                                                                                                               dist_x,
                                                                                                                               dist_y,
                                                                                                                               dist_z)
                
                # Getting the left (lower) triangle of H and S matrix
                self.H_dftb[jtotal_orb][itotal_orb] = self.H_dftb[itotal_orb][jtotal_orb]
                self.S_dftb[jtotal_orb][itotal_orb] = self.S_dftb[itotal_orb][jtotal_orb]

    def print_H_S_matrix(self):
                          
        # Printing the Hamiltonian matrix
        fileOut = file("dftb.out", "w")
        fileOut.write("DFTB matrix:")
        fileOut.write("\n")
        fileOut.write("------------" + "\n")
   
        count = 0
        for i in xrange(self.total_orb):
            if count < self.total_orb:
               fileOut.write("%9s" %  (i + 1) + "  ")
            else:
               fileOut.write("%s" % (i + 1) + "  " + "\n")
               count = 0
        fileOut.write("\n")

        count = 0
        index = 0
        for iH_dftb in self.H_dftb:
            index += 1
            fileOut.write("%s %s %s" % (index, self.atom_symb[index - 1],  self.atom_orb_sym[index - 1]) + "  ")
            for element in iH_dftb:
                count += 1
                if count < len(iH_dftb):
                   fileOut.write("%f" % (element) + "  ",)
                else:
                  fileOut.write("%f" % (element) + "  " + "\n")
                  count = 0

        # Printing the Overlap matrix
        fileOut.write("\n")
        fileOut.write("Overlap matrix:" + "\n")
        fileOut.write("---------------" + "\n")

        count = 0
        for i in xrange(self.total_orb):
            if count < self.total_orb:
               fileOut.write("%9s" %  (i + 1) + "  ")
            else:
               fileOut.write("%s" % (i + 1) + "  " + "\n")
               count = 0
        fileOut.write("\n")

        count = 0
        index = 0
        for iS_dftb in self.S_dftb:
            index += 1
            fileOut.write("%s %s %s" % (index, self.atom_symb[index - 1],  self.atom_orb_sym[index - 1]) + "  ")
            for element in iS_dftb:
                count += 1
                if count < len(iS_dftb):
                   fileOut.write("%f" % (element) + "  ",)
                else:
                  fileOut.write("%f" % (element) + "  " + "\n")
                  count = 0

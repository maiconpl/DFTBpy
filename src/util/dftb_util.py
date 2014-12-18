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

import sys
from math import sqrt
from model.fit import Fit

class DFTBMatrixUtil(object):

    def __init__(self, slako_read_list, 
                 atom_type, 
                 n_atom,
                 total_charge, 
                 max_angular_momentum, 
                 max_all_momentum, 
                 xyz_atom_symbols, 
                 xyz,
                 debug):

        self.slako_read_list = slako_read_list
        self.atom_type = atom_type
        self.n_atom = n_atom
        self.total_charge = total_charge
        self.max_angular_momentum = max_angular_momentum
        self.max_all_momentum = max_all_momentum
        self.xyz_atom_symbols = xyz_atom_symbols
        self.xyz = xyz
        self.debug = debug

        self.atom_orb_sym= []  # list of the symmetry of the atomic orbitals.
        self.H_diag_dftb = []  # list of the diagonal H.
        self.S_diag_dftb = []  # list of the diagonal S.
        self.total_orb = 0
        self.atom_symb = []  # symbol of the atoms related with the angular momentum list.
        self.atom_numb = []  # number of the atom related with the angular momentum list.
        self.atom_type_numb = []  # list related with the atoms type and the angular momentum list.
        self.atoms_type_combinations_list = [] # list of combinations of atoms type symbols
        self.dist_matrix = []
        self.dist_x = []
        self.dist_y = []
        self.dist_z = []
        self.orbitals_index_range_list_per_atom = []
        
        self.fit = Fit()
                        
    def build_data_lists(self):

        # This routine coorrelate the symmetry of the Hamiltonian and Overlap integrals.
        atom_orb_label = []  # atoms and orbitals label.
        count = 0
        tmp = []
        n_atom_type = len(self.atom_type)
                
        for iatom in xrange(self.n_atom):
            for jatom_type in xrange(n_atom_type):

                # creating the variable "slako_read" related with the previous object: slako_read_list
                slako_read = self.slako_read_list[jatom_type][jatom_type]

                if self.xyz_atom_symbols[iatom] == self.atom_type[jatom_type] and self.max_angular_momentum[jatom_type] == "s":
                   count = count +  1
                   # creating a list with the symbol of the valence orbitals:
                   tmp = "s"
      	           self.atom_orb_sym.append(tmp)
                   # colect the eigenvalues of s orbital:
                   tmp2 = slako_read.orb_energy[2]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # creating a list with the symbol of the atoms related with the angular momentum s:
                   tmp4 = self.xyz_atom_symbols[iatom]
                   self.atom_symb.append(tmp4)
                   # creating a list with the number of the atoms related with the angular momentum s:
                   tmp5 = iatom
                   self.atom_numb.append(tmp5)
                   # creating a list with the number of the atoms related with the angular momentum s:
                   tmp6 = jatom_type
                   self.atom_type_numb.append(tmp6)
                   #
                if self.xyz_atom_symbols[iatom] == self.atom_type[jatom_type] and self.max_angular_momentum[jatom_type] == "p":
                   count = count +  4
                   # creating a list with the symbol of the valence orbitals:
                   tmp = "s"
      	           self.atom_orb_sym.append(tmp)
                   #
                   tmp = "px"
        	   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "py"
         	   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "pz"
      	           self.atom_orb_sym.append(tmp)
                   # colect the eigenvalues of s orbital:
                   tmp2 = slako_read.orb_energy[2]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # colect the eigenvalues of p orbital:
                   tmp2 = slako_read.orb_energy[1]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # creating a list with the symbol of the atoms related with the angular momentum s and p:
                   tmp4 = self.xyz_atom_symbols[iatom]
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   # creating a list with the number of the atoms related with the angular momentum s and p:
                   tmp5 = iatom
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   # creating a list with the number of the atoms related with the angular momentum s and p:
                   tmp6 = jatom_type
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   #
                if self.xyz_atom_symbols[iatom] == self.atom_type[jatom_type] and self.max_angular_momentum[jatom_type] == "d":
                   count = count +  9
                   # creating a list with the symbol of the valence orbitals:
                   tmp = "s"
                   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "px"
                   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "py"
                   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "pz"
                   self.atom_orb_sym.append(tmp)
                   #
                   tmp = "dxy"
      	           self.atom_orb_sym.append(tmp)
                   #
                   tmp = "dyz"
      	           self.atom_orb_sym.append(tmp)
                   #
                   tmp = "dxz"
      	           self.atom_orb_sym.append(tmp)
                   #
                   tmp = "dx2-y2"
      	           self.atom_orb_sym.append(tmp)
                   #
                   tmp = "dz2"
      	           self.atom_orb_sym.append(tmp)
                   # colect the eigenvalues of s orbital:
                   tmp2 = slako_read.orb_energy[2]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # colect the eigenvalues of p orbital:
                   tmp2 = slako_read.orb_energy[1]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # colect the eigenvalues of d orbital:
                   tmp2 = slako_read.orb_energy[0]
                   tmp3 = 1.0
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   #
                   self.H_diag_dftb.append(tmp2)
                   self.S_diag_dftb.append(tmp3)
                   # creating a list with the symbol of the atoms related with the angular momentum s and p and d:
                   tmp4 = self.xyz_atom_symbols[iatom]
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   self.atom_symb.append(tmp4)
                   # creating a list with the number of the atoms related with the angular momentum s and p and d:
                   tmp5 = iatom
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   self.atom_numb.append(tmp5)
                   # creating a list with the number of the atoms related with the angular momentum s p and d:
                   tmp6 = jatom_type
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   self.atom_type_numb.append(tmp6)
                   #
        self.total_orb = self.total_orb + count
        count = 0
# Oy for debugging purposes
#       print "Diagonal elements:"
#       print self.H_diag_dftb
#       print "Symmetry of the valence atomic orbitals:"
#       print self.atom_orb_sym
#       print "Symbol of the atoms:"
#        print self.atom_symb
#       print "Symbol of the number of the atoms which apear, respectively, in the xyz file:"
#        print self.atom_numb
#       print "List with the number of the atoms related with the atoms type."
#        print  self.atom_type_numb
#       print "Total number of valence orbitals:", self.total_orb

        #  Criating a list with the total number of valence orbitals
#       for i in xrange(self.total_orb):
#           tmp = i
#           atom_orb_label.append(tmp)
#       print atom_orb_label

      #  Getting the total number of electrons:

    def build_orbitals_index_range_list_per_atom(self):
    
        initial_index = 0
        final_index = 0

        for iAtom in xrange(self.n_atom):

            count = 0
            isIndex = False
            for i in xrange(self.total_orb):
                if self.atom_numb[i] == iAtom:
                   count = count + 1
                   if isIndex == False:
                      initial_index = i
                      isIndex = True

            final_index = initial_index + count - 1

            self.orbitals_index_range_list_per_atom.append((initial_index, final_index))

        return self.orbitals_index_range_list_per_atom
    
    def build_H_S_dist_list(self):
        
        # Defining the valid region of interpolation
        distance_grid = []
        integral_H_grid = []
        integral_S_grid = []
        line_dist = 0 # relationship between the distance and the number of the line.
        isTrue = False
        
        n_atom_type = len(self.atom_type)
        H_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]
        S_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]
        dist_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]
        dist_step_list = [[0.0 for i in xrange(n_atom_type)] for j in xrange(n_atom_type)]

        for iatom_type in xrange(n_atom_type):
            for jatom_type in xrange(n_atom_type):

                # creating the variable "slako_read" related with the previous object: slako_read_list
                slako_read = self.slako_read_list[iatom_type][jatom_type]

                dist_step_list[iatom_type][jatom_type] = slako_read.dist_step 

                for isymmetry in xrange(10):

                    line_dist = 0 # relationship between the distance and the number of the line.
#                    isTrue = False
                    
                    distance_grid = []
                    integral_H_grid = []
                    integral_S_grid = []
                    
                    for igrid in xrange(len(slako_read.skf_H)):
                        
                        line_dist = line_dist + 1

#                         if float(slako_read.skf_H[igrid][isymmetry]) == 0.0 and isTrue == False:
#                             pass
#                         if float(slako_read.skf_H[igrid][isymmetry]) != 0.0:
#                             isTrue = True
                        distance_grid.append(float(slako_read.dist_step*line_dist))
                        integral_H_grid.append(float(slako_read.skf_H[igrid][isymmetry]))
                        integral_S_grid.append(float(slako_read.skf_S[igrid][isymmetry]))

                    dist_list[iatom_type][jatom_type][isymmetry] = distance_grid    
                    H_list[iatom_type][jatom_type][isymmetry] = integral_H_grid
                    S_list[iatom_type][jatom_type][isymmetry] = integral_S_grid
#                     print iatom_type, jatom_type, isymmetry, len(distance_grid), len(integral_H_grid), len(integral_S_grid)
#                     print "za"
#                     print H_list[iatom_type][jatom_type][isymmetry], len(H_list[iatom_type][jatom_type][isymmetry])
#                     print S_list[iatom_type][jatom_type][isymmetry], len(S_list[iatom_type][jatom_type][isymmetry])
#                     print dist_list[iatom_type][jatom_type][isymmetry], len(dist_list[iatom_type][jatom_type][isymmetry])

        return H_list, S_list, dist_list, dist_step_list

    def build_skf_fit_list(self, dist_list, H_list, S_list):

        n_atom_type = len(self.atom_type)
        fit_H_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]
        fit_S_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]

        for iatom_type in xrange(n_atom_type):
            for jatom_type in xrange(n_atom_type):

                for isymmetry in xrange(10):

                    distance_grid = dist_list[iatom_type][jatom_type][isymmetry]
                    fit_H_list[iatom_type][jatom_type][isymmetry] = self.fit.slako_integrals_fit2(distance_grid, 
                                                                                             H_list[iatom_type][jatom_type][isymmetry])
                    
                    fit_S_list[iatom_type][jatom_type][isymmetry] = self.fit.slako_integrals_fit2(distance_grid,
                                                                                             S_list[iatom_type][jatom_type][isymmetry])

        return fit_H_list, fit_S_list

    def print_skf_fit(self, dist_list, H_list, S_list, fit_H_list, fit_S_list):

        import matplotlib.pyplot as plt
        import numpy as np
        
        symmetry = 7
        
        x = dist_list[0][0][symmetry]
        y = H_list[0][0][symmetry]
        xnew = xnew = np.arange(0.5, len(x), 0.2321)
        ynew_spline = fit_H_list[0][0][symmetry](xnew)
        
#         n_atom_type = len(self.atom_type)
#         for iatom_type in xrange(n_atom_type):
#             for jatom_type in xrange(n_atom_type):
#         
#                 fit_H_list[iatom_type][jatom_type][9], dist_list[iatom_type][jatom_type][isymmetry]
                
                
        plt.figure()
#        plt.plot(x, y, 'x', xnew, ynew_spline, '*', xnew, np.sin(xnew), xnew, ynew_interp1d, 'o')
        plt.plot(x, y, '-', xnew, ynew_spline, '*')#, xnew, np.sin(xnew), xnew, ynew_interp1d, 'o')

        plt.legend(['Real', 'Interpolated'])
        plt.axis([dist_list[0][0][9][0] - 1.0, dist_list[0][0][9][-1] + 0.4, H_list[0][0][9][0] - 0.1, H_list[0][0][9][-1] + 1.0])
        plt.title('InterpolatedUnivariateSpline <s|H|s>')
        plt.show()

    def test_polint_fit_(self, dist_list, dist_step, H_list, S_list):
                
        dist_atoms = 19.642812341234
        
        a = self.fit.slako_integrals_fit_from_C(dist_list, dist_step, H_list, dist_atoms)
        
#        print dist_atoms, a 

    def build_integral_cutoff_list(self, H_list, dist_list):

        # We evaluate the matrix elements cutoff just for the element matrix: <s|H|s>.
        
        n_atom_type = len(self.atom_type)
        H_cutoff_list = [[ [0.0 for i in xrange(10)] for j in xrange(n_atom_type)] for k in xrange(n_atom_type)]
        inferior_distance_cutoff = 0.0

        for iatom_type in xrange(n_atom_type):
            for jatom_type in xrange(n_atom_type):

                for isymmetry in xrange(10):

                    for i in xrange(len(dist_list[iatom_type][jatom_type][isymmetry])):
                        if H_list[iatom_type][jatom_type][isymmetry][i] != 0.0:
                            inferior_distance_cutoff = dist_list[iatom_type][jatom_type][isymmetry][i]
                            break

                    superior_distance_cutoff = dist_list[iatom_type][jatom_type][isymmetry][-1]
                    H_cutoff_list[iatom_type][jatom_type][isymmetry] = inferior_distance_cutoff, superior_distance_cutoff
                    
        return H_cutoff_list

    def get_number_of_electrons(self):

        self.n_electrons = [] # total number of electrons
        self.n_electrons_atom = [] # number of valence electrons per atom of the system
        n_atom_type = n_atom_type = len(self.atom_type)
        count = 0
        for iatom in xrange(self.n_atom):
            for jatom_type in xrange(n_atom_type):
                 
                # creating the variable "slako_read" related with the previous object: slako_read_list
                slako_read = self.slako_read_list[jatom_type][jatom_type]
                #
                if self.xyz_atom_symbols[iatom] == self.atom_type[jatom_type]:
                   count += 1
                   self.n_electrons_atom.append(sum(slako_read.orb_ocupation))

        self.n_electrons = int(sum(self.n_electrons_atom) - self.total_charge)
        
        if self.n_electrons < 0:
            
            print ""
            print "ERROR!!!"
            print "The total number of electrons " + "'" + str(self.n_electrons) + "' is negative."
            print "The charge " + "'" + str(self.total_charge) + "' " + "can not be greater than the number of valence electrons"
            print "of the atom(s) of the system. In this case " + "'" + str(sum(self.n_electrons_atom)) + "' electrons."
            sys.exit()

        if (self.n_electrons % 2) != 0:

            print ""
            print "ERROR!!!"
            print "The total number of electrons is odd: " + str(self.n_electrons) + "."
            print "It should be even. We do not deal with DFTB/SCC-DFTB open shell calculations yet."
            print "Please, check the atom(s) in the systems and/or the set charge in the input."
            sys.exit()

    def build_dist_matrix(self, n_atom, xyz):

        dist_matrix = [ [0.0 for i in xrange(n_atom)] for j in xrange(n_atom)] # defining the dimension of the matrix.
        dist_x = [ [0.0 for i in xrange(n_atom)] for j in xrange(n_atom)]
        dist_y = [ [0.0 for i in xrange(n_atom)] for j in xrange(n_atom)]
        dist_z = [ [0.0 for i in xrange(n_atom)] for j in xrange(n_atom)]

        for iatom in xrange(n_atom):
            for jatom in xrange(n_atom):
                dist_x[iatom][jatom] = (xyz[iatom][0] - xyz[jatom][0])
                dist_y[iatom][jatom] = (xyz[iatom][1] - xyz[jatom][1])
                dist_z[iatom][jatom] = (xyz[iatom][2] - xyz[jatom][2])
                dist_matrix[iatom][jatom] = sqrt( (dist_x[iatom][jatom])**2 + (dist_y[iatom][jatom])**2 + (dist_z[iatom][jatom])**2)
                dist_matrix[jatom][iatom] = dist_matrix[iatom][jatom]
                
        return dist_matrix, dist_x, dist_y, dist_z

    def print_dist_matrix(self):

        print "----------------------------"
        print "Interatomic distances matrix"
        print "----------------------------"
     
        count = 0
        for idist_matrix in self.dist_matrix:
            for element in idist_matrix:
                count += 1
                if count < len(idist_matrix):
                   print "%f" % (element) + "  ",
                else:
                   print "%f" % (element) + ""
                   count = 0

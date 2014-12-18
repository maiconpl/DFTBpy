# Created by MPL on Oct 24, 2014. Last modification
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

import os
from math import sqrt
from model.slater_koster_transformations import SlaterKosterTransformations
_path = os.path.dirname(os.path.realpath(__file__)) + "/" + "../lib/build_dgamma_matrix.so"

from ctypes import cdll, byref, c_int, c_double, c_char
lib_dgamma_matrix = cdll.LoadLibrary(_path)

class ForcesOnAtom(object):

    def __init__(self, dist_list,
                 dist_step_list,
                 H_list,
                 S_list,
                 slako_read_list,
                 dist_x,
                 dist_y,
                 dist_z,
                 dist_matrix,
                 n_atom,
                 xyz_atom_symbols,
                 atom_type,
                 atom_symb,
                 atom_numb,
                 atom_type_numb,
                 atoms_type_combinations_list,
                 orbitals_index_range_list_per_atom,
                 max_all_momentum,
                 atom_orb_sym,
                 total_orb,
                 P_dftb,
                 P_dftb_e,
                 net_charge,
                 scc_matrix_without_S,
                 isSCC,
                 debug):

        self.dist_list = dist_list
        self.dist_step_list = dist_step_list
        self.H_list = H_list
        self.S_list = S_list
        self.slako_read_list = slako_read_list
        self.dist_x = dist_x
        self.dist_y = dist_y
        self.dist_z = dist_z
        self.dist_matrix = dist_matrix
        self.n_atom = n_atom
        self.xyz_atom_symbols = xyz_atom_symbols
        self.atom_type = atom_type
        self.atom_symb = atom_symb
        self.atom_numb = atom_numb
        self.atoms_type_combinations_list = atoms_type_combinations_list
        self.orbitals_index_range_list_per_atom = orbitals_index_range_list_per_atom
        self.atom_type_numb = atom_type_numb
        self.max_all_momentum = max_all_momentum
        self.atom_orb_sym = atom_orb_sym
        self.total_orb = total_orb
        self.P_dftb = P_dftb
        self.P_dftb_e = P_dftb_e
        self.net_charge = net_charge
        self.isSCC = isSCC
        self.debug = debug
        
        self.hubbard_s = []
        self.dgamma = []

        self.H_plus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.H_plus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.H_plus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        self.S_plus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.S_plus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.S_plus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        self.H_minus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.H_minus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.H_minus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        
        self.S_minus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.S_minus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.S_minus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        self.electronic_force_x = [0.0 for i in xrange(self.n_atom)]
        self.electronic_force_y = [0.0 for i in xrange(self.n_atom)]
        self.electronic_force_z = [0.0 for i in xrange(self.n_atom)]

        self.repulsion_force_x = [0.0 for i in xrange(self.n_atom)]
        self.repulsion_force_y = [0.0 for i in xrange(self.n_atom)]
        self.repulsion_force_z = [0.0 for i in xrange(self.n_atom)]

        self.scc_force_x = [0.0 for i in xrange(self.n_atom)]
        self.scc_force_y = [0.0 for i in xrange(self.n_atom)]
        self.scc_force_z = [0.0 for i in xrange(self.n_atom)]

        self.sk_transformations = SlaterKosterTransformations(dist_list,
                                                              dist_step_list,
                                                              H_list,
                                                              S_list,
                                                              max_all_momentum)
        
        if self.isSCC == True:

            self.scc_matrix_without_S = scc_matrix_without_S
        
            # Defining the ctypes

            self._n_atom = c_int(self.n_atom)
            self._total_orb = c_int(total_orb)
            _dgamma = (c_double*(self.n_atom))*self.n_atom
            self._dgamma = _dgamma()
            _scc_matrix = (c_double*(total_orb))*total_orb
            self._scc_matrix = _scc_matrix()

            _dist_x = (c_double*(self.n_atom))*self.n_atom
            self._dist_x = _dist_x()
            _dist_y = (c_double*(self.n_atom))*self.n_atom
            self._dist_y = _dist_y()
            _dist_z = (c_double*(self.n_atom))*self.n_atom
            self._dist_z = _dist_z()
            _dist_matrix = (c_double*(self.n_atom))*self.n_atom
            self._dist_matrix = _dist_matrix()
            _xyz_atom_symbols = ((c_char*3)*self.n_atom)
            self._xyz_atom_symbols = _xyz_atom_symbols()
            _hubbard_s = c_double*(self.n_atom)
            self._hubbard_s = _hubbard_s()
    
            for i in xrange(self.n_atom):
                self._hubbard_s[i] = 0.0
                self._xyz_atom_symbols[i].value = xyz_atom_symbols[i]
                for j in xrange(i, self.n_atom):
                    self._dgamma[i][j] = 0.0
                    self._dgamma[j][i] = 0.0
                    self._dist_x[i][j] = dist_x[i][j]
                    self._dist_x[j][i] = dist_x[i][j]
                    self._dist_y[i][j] = dist_y[i][j]
                    self._dist_y[j][i] = dist_y[i][j]
                    self._dist_z[i][j] = dist_z[i][j]
                    self._dist_z[j][i] = dist_z[i][j]
                    self._dist_matrix[i][j] = dist_matrix[i][j]
                    self._dist_matrix[j][i] = dist_matrix[i][j]
            
    def get_electronic_forces(self):
 
#         print "xyz:"
#         print self.xyz
#         
#         print "xyz_atom_symbols"
#         print self.xyz_atom_symbols
#         
#         print "atom_symb"
#         print self.atom_symb
# 
#         print "atom_numb"
#         print self.atom_numb
#         
#         print "atom_type_numb"
#         print self.atom_type_numb
#         
#         print "max_all_momentum"
#         print self.max_all_momentum
#         
#         print "atom_orb_sym"
#         print self.atom_orb_sym
#         
#         print "total_orb"
#         print self.total_orb
# #        print self.P_dftb

#        step = 0.00188972599 # step in bohr: 0.1 A = 0.188972599 bohr
        step = 0.0
        der_H_x = 0.0
        der_H_y = 0.0
        der_H_z = 0.0
        der_S_x = 0.0
        der_S_y = 0.0
        der_S_z = 0.0
        
        isDerivativeDone = [ [False for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        dist_matrix_plus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        dist_matrix_plus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        dist_matrix_plus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        dist_matrix_minus_x = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        dist_matrix_minus_y = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        dist_matrix_minus_z = [ [0.0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        
        sequence_list = list(range(0, self.total_orb, 1))

        for iAtom in xrange(self.n_atom):

#             print "Analysis:", self.xyz_atom_symbols[iAtom]

            der_H_x = 0.0
            der_H_y = 0.0
            der_H_z = 0.0
            der_S_x = 0.0
            der_S_y = 0.0
            der_S_z = 0.0

	    # This is a highly optimized set of nested loops.
            # But, it makes the code less clear.
            # All redundancies, including an important if, were removed just by exploring the lists
            # 'self.orbitals_index_range_list_per_atom' and 'sequence_list'
            # to mark the index to be evaluated.
            # The copy 'forces.py.cp31' is a backup with a more clear code.
            for j in xrange(self.orbitals_index_range_list_per_atom[iAtom][0],
                            self.orbitals_index_range_list_per_atom[iAtom][1] + 1):

                k = 0
                isTrue = False
                size_to_not_consider = (self.orbitals_index_range_list_per_atom[iAtom][1] + 1) - self.orbitals_index_range_list_per_atom[iAtom][0]
                for l in xrange(self.total_orb - size_to_not_consider):
                    
                    if l < self.orbitals_index_range_list_per_atom[iAtom][0]:
                        k = sequence_list[l]
                    
                    elif isTrue == False:
                        k = sequence_list[self.orbitals_index_range_list_per_atom[iAtom][1] + 1]

                        isTrue = True
                    else:
                        k = k + 1

#                    print "iAtom", iAtom, j, k, self.atom_numb[j], self.atom_numb[k]
                    
#                    if j != k and self.dist_matrix[self.atom_numb[j]][self.atom_numb[k]] != 0.0 and iAtom == self.atom_numb[j]:

#                    print "iAtom2", iAtom, j, k, self.atom_numb[j], self.atom_numb[k]
                    #print "->", self.atom_symb[j], self.atom_symb[k], self.atom_numb[j], self.atom_numb[k], self.dist_matrix[j][k], self.atom_orb_sym[j],  self.atom_orb_sym[k], self.xyz[self.atom_numb[j]][0], self.xyz[self.atom_numb[j]][1], self.xyz[self.atom_numb[j]][2]

                    step = 0.1*self.slako_read_list[self.atom_type_numb[j]][self.atom_type_numb[k]].dist_step

                    if isDerivativeDone[j][k] == False: # This simple 'if' is for optimization purposes, it yields ~33% of performance.

                        isDerivativeDone[j][k] = True
                        isDerivativeDone[k][j] = True

                        dist_x = self.dist_x[self.atom_numb[j]][self.atom_numb[k]]
                        dist_y = self.dist_y[self.atom_numb[j]][self.atom_numb[k]]
                        dist_z = self.dist_z[self.atom_numb[j]][self.atom_numb[k]]
        
                        integral_symmetry = self.atom_orb_sym[j] + self.atom_orb_sym[k]
                        
                        dist_x_plus = dist_x + step
                        dist_y_plus = dist_y + step
                        dist_z_plus = dist_z + step
            
                        dist_x_minus = dist_x - step
                        dist_y_minus = dist_y - step
                        dist_z_minus = dist_z - step
        
    #                    H, S = self.sk_transformations.get_skf(self.atom_type_numb[j],
    #                                                           self.atom_type_numb[k],
    #                                                           integral_symmetry,
    #                                                           self.dist_matrix[self.atom_numb[j]][self.atom_numb[k]],
    #                                                           dist_x,
    #                                                           dist_y,
    #                                                           dist_z)
        
                        # H_plus_x, S_plus_x
                        dist_matrix_plus_x[j][k] = sqrt(dist_x_plus**2.0 + dist_y**2.0 + dist_z**2.0)
                        dist_atoms = dist_matrix_plus_x[j][k]
                        self.H_plus_x[j][k], self.S_plus_x[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                             self.atom_type_numb[k],
                                                                             integral_symmetry,
                                                                             dist_atoms,
                                                                             dist_x_plus,
                                                                             dist_y,
                                                                             dist_z)
                        self.H_plus_x[k][j] = -self.H_plus_x[j][k]
                        self.S_plus_x[k][j] = -self.S_plus_x[j][k]
        
                        # H_plus_y, S_plus_y
                        dist_matrix_plus_y[j][k] = sqrt(dist_x**2.0 + dist_y_plus**2.0 + dist_z**2.0)
                        dist_atoms = dist_matrix_plus_y[j][k]
                        self.H_plus_y[j][k], self.S_plus_y[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                             self.atom_type_numb[k],
                                                                             integral_symmetry,
                                                                             dist_atoms,
                                                                             dist_x,
                                                                             dist_y_plus,
                                                                             dist_z)
                        self.H_plus_y[k][j] = -self.H_plus_y[j][k]
                        self.S_plus_y[k][j] = -self.S_plus_y[j][k]
        
                        # H_plus_z, S_plus_z
                        dist_matrix_plus_z[j][k] = sqrt(dist_x**2.0 + dist_y**2.0 + dist_z_plus**2.0)
                        dist_atoms = dist_matrix_plus_z[j][k]
                        self.H_plus_z[j][k], self.S_plus_z[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                             self.atom_type_numb[k],
                                                                             integral_symmetry,
                                                                             dist_atoms,
                                                                             dist_x,
                                                                             dist_y,
                                                                             dist_z_plus)
                        self.H_plus_z[k][j] = -self.H_plus_z[j][k]
                        self.S_plus_z[k][j] = -self.S_plus_z[j][k]
        
                        # H_minus_x, S_minus_x
                        dist_matrix_minus_x[j][k] = sqrt(dist_x_minus**2.0 + dist_y**2.0 + dist_z**2.0)
                        dist_atoms = dist_matrix_minus_x[j][k]
                        self.H_minus_x[j][k], self.S_minus_x[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                               self.atom_type_numb[k],
                                                                               integral_symmetry,
                                                                               dist_atoms,
                                                                               dist_x_minus,
                                                                               dist_y,
                                                                               dist_z)
                        self.H_minus_x[k][j] = -self.H_minus_x[j][k]
                        self.S_minus_x[k][j] = -self.S_minus_x[j][k]
                        
                        # H_minus_y, S_minus_y
                        dist_matrix_minus_y[j][k] = sqrt(dist_x**2.0 + dist_y_minus**2.0 + dist_z**2.0)
                        dist_atoms = dist_matrix_minus_y[j][k]
                        self.H_minus_y[j][k], self.S_minus_y[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                               self.atom_type_numb[k],
                                                                               integral_symmetry,
                                                                               dist_atoms,
                                                                               dist_x,
                                                                               dist_y_minus,
                                                                               dist_z)
                        self.H_minus_y[k][j] = -self.H_minus_y[j][k]
                        self.S_minus_y[k][j] = -self.S_minus_y[j][k]
                        
                        # H_minus_z, S_minus_z
                        dist_matrix_minus_z[j][k] = sqrt(dist_x**2.0 + dist_y**2.0 + dist_z_minus**2.0)
                        dist_atoms = dist_matrix_minus_z[j][k]
                        self.H_minus_z[j][k], self.S_minus_z[j][k] = self.sk_transformations.get_skf(self.atom_type_numb[j],
                                                                               self.atom_type_numb[k],
                                                                               integral_symmetry,
                                                                               dist_atoms,
                                                                               dist_x,
                                                                               dist_y,
                                                                               dist_z_minus)
                        self.H_minus_z[k][j] = -self.H_minus_z[j][k]
                        self.S_minus_z[k][j] = -self.S_minus_z[j][k]

#                    print "H, S"
#                    print  H, S
#                    print "H_plus_x, H_minus_x, S_plus_x, S_minus_x"
#                    print  H_plus_x, H_minus_x, S_plus_x, S_minus_x
#                    print "H_plus_y, H_minus_y, S_plus_y, S_minus_y"
#                    print  H_plus_y, H_minus_y, S_plus_y, S_minus_y
#                    print "H_plus_z, H_minus_z, S_plus_z, S_minus_z"
#                    print  H_plus_z, H_minus_z, S_plus_z, S_minus_z

                    # Puv = 2.0*sum_i^(N/2)*CuiCvi
                    der_H_x = der_H_x - self.P_dftb[j][k]*((self.H_plus_x[j][k] - self.H_minus_x[j][k])/(2.0*step))
                    der_H_y = der_H_y - self.P_dftb[j][k]*((self.H_plus_y[j][k] - self.H_minus_y[j][k])/(2.0*step))
                    der_H_z = der_H_z - self.P_dftb[j][k]*((self.H_plus_z[j][k] - self.H_minus_z[j][k])/(2.0*step))
    
                    der_S_x = der_S_x + self.P_dftb_e[j][k]*((self.S_plus_x[j][k] - self.S_minus_x[j][k])/(2.0*step))
                    der_S_y = der_S_y + self.P_dftb_e[j][k]*((self.S_plus_y[j][k] - self.S_minus_y[j][k])/(2.0*step))
                    der_S_z = der_S_z + self.P_dftb_e[j][k]*((self.S_plus_z[j][k] - self.S_minus_z[j][k])/(2.0*step))
                    
                    if self.isSCC == True:

                        der_S_x = der_S_x - self.P_dftb[j][k]*self.scc_matrix_without_S[j][k]*((self.S_plus_x[j][k] - self.S_minus_x[j][k])/(2.0*step))
                        der_S_y = der_S_y - self.P_dftb[j][k]*self.scc_matrix_without_S[j][k]*((self.S_plus_y[j][k] - self.S_minus_y[j][k])/(2.0*step))
                        der_S_z = der_S_z - self.P_dftb[j][k]*self.scc_matrix_without_S[j][k]*((self.S_plus_z[j][k] - self.S_minus_z[j][k])/(2.0*step))
    
            self.electronic_force_x[iAtom] = -2.0*(der_H_x + 1.0*der_S_x)
            self.electronic_force_y[iAtom] = -2.0*(der_H_y + 1.0*der_S_y)
            self.electronic_force_z[iAtom] = -2.0*(der_H_z + 1.0*der_S_z)
            # The multiplication by 2.0 is because we worked with just half (1/2) part of the
            # H and S matrix to evaluate the forces on atoms due to the if condition.

        if self.debug == True:
            print ""
            print "Electronic forces:"
            for iAtom in xrange(self.n_atom):
                print ("%2s % .12f % .12f % .12f, module:% .12f" % (self.xyz_atom_symbols[iAtom], \
                                                                    self.electronic_force_x[iAtom], \
                                                                    self.electronic_force_y[iAtom], \
                                                                    self.electronic_force_z[iAtom], \
                                                                    sqrt(self.electronic_force_x[iAtom]**2 + self.electronic_force_y[iAtom]**2 + self.electronic_force_z[iAtom]**2)))
    
    def get_repulsion_forces(self):

        n_atom_type = len(self.atom_type)
        
        force_x = 0.0
        force_y = 0.0
        force_z = 0.0

        for iatom in xrange(self.n_atom):

            force_x = 0.0
            force_y = 0.0
            force_z = 0.0

            for jatom in xrange(self.n_atom):
                
                poly_rep_force = 0.0

                if iatom != jatom:
                    
                    for i in xrange(n_atom_type):
                        for j in xrange(n_atom_type):
                            if self.xyz_atom_symbols[iatom] + self.xyz_atom_symbols[jatom] == \
                                self.atoms_type_combinations_list[i][j]:
                                slako_read = self.slako_read_list[i][j]
                                break
    
                    atoms_dist = 0.0
                    atoms_dist = self.dist_matrix[iatom][jatom]
    
                    diff = slako_read.cut_radius - atoms_dist
    
                    if atoms_dist <= slako_read.cut_radius:
                        poly_rep_force = 2.0*slako_read.poly_coef[0]*diff + 3.0*slako_read.poly_coef[1]*(diff**2) + \
                                   4.0*slako_read.poly_coef[2]*(diff**3) + 5.0*slako_read.poly_coef[3]*(diff**4) + \
                                   6.0*slako_read.poly_coef[4]*(diff**5) + 7.0*slako_read.poly_coef[5]*(diff**6) + \
                                   8.0*slako_read.poly_coef[6]*(diff**7)
                    
                    else:
                        poly_rep_force = 0.0
    
                    force_x = force_x - (poly_rep_force)*(self.dist_x[iatom][jatom]/atoms_dist)
                    force_y = force_y - (poly_rep_force)*(self.dist_y[iatom][jatom]/atoms_dist)
                    force_z = force_z - (poly_rep_force)*(self.dist_z[iatom][jatom]/atoms_dist)

            self.repulsion_force_x[iatom] = force_x
            self.repulsion_force_y[iatom] = force_y
            self.repulsion_force_z[iatom] = force_z

        if self.debug == True:                
            print "Repulsion forces:"
            for iatom in xrange(self.n_atom):
                print ("%2s % .12f % .12f % .12f, module:% .12f" % (self.xyz_atom_symbols[iatom],  \
                                                                    self.repulsion_force_x[iatom], \
                                                                    self.repulsion_force_y[iatom], \
                                                                    self.repulsion_force_z[iatom], 
                                                                    sqrt(self.repulsion_force_x[iatom]**2 + self.repulsion_force_y[iatom]**2 + self.repulsion_force_z[iatom]**2)))

    def get_hubbard_parameter(self):

        n_atom_type = len(self.atom_type)
        self.hubbard_s = [0 for i in xrange(self.n_atom)]

        for iatom in xrange(self.n_atom):

            for i in xrange(n_atom_type):
                if self.xyz_atom_symbols[iatom] + self.xyz_atom_symbols[iatom] == \
                    self.atoms_type_combinations_list[i][i]:
                    slako_read = self.slako_read_list[i][i]
                    break

            self.hubbard_s[iatom] = slako_read.hubbard_parameter[2]

    def get_gamma_forces(self):

        self.get_hubbard_parameter()

        for i in xrange(self.n_atom):
            self._hubbard_s[i] = self.hubbard_s[i]

        self.dgamma = [ [0.0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.

        lib_dgamma_matrix.build_dgamma_matrix(self._n_atom,
                                   byref(self._dist_matrix),
                                   byref(self._xyz_atom_symbols),
                                   byref(self._hubbard_s),
                                   byref(self._dgamma))

        self.dgamma = [ [0.0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.

        for i in xrange(self.n_atom):
            for j in xrange(i, self.n_atom):
                self.dgamma[i][j] = self._dgamma[i][j]
                self.dgamma[j][i] = self._dgamma[i][j]

    def get_scc_forces(self):

        self.get_gamma_forces()

        force_x = 0.0
        force_y = 0.0
        force_z = 0.0

        for i in xrange(self.n_atom):

            force_x = 0.0
            force_y = 0.0
            force_z = 0.0

            for j in xrange(self.n_atom):

                dist_x = self.dist_x[i][j]
                dist_y = self.dist_y[i][j]
                dist_z = self.dist_z[i][j]

                if i != j:
            
                   tmp_scc = 0.0
                   tmp_scc = -self.dgamma[i][j]*(-self.net_charge[j])

                   force_x = force_x - (tmp_scc*dist_x)
                   force_y = force_y - (tmp_scc*dist_y)
                   force_z = force_z - (tmp_scc*dist_z)

            self.scc_force_x[i] = force_x*(-self.net_charge[i])
            self.scc_force_y[i] = force_y*(-self.net_charge[i])
            self.scc_force_z[i] = force_z*(-self.net_charge[i])

        if self.debug == True:
            print "SCC forces:"
            for iatom in xrange(self.n_atom):
                print ("%2s % .12f % .12f % .12f, module:% .12f" % (self.xyz_atom_symbols[iatom], \
                                                                    self.scc_force_x[iatom], \
                                                                    self.scc_force_y[iatom], \
                                                                    self.scc_force_z[iatom], \
                                                                    sqrt(self.scc_force_x[iatom]**2 + self.scc_force_y[iatom]**2 + self.scc_force_z[iatom]**2)))

    def get_dftb_forces(self):

        self.get_electronic_forces()
        self.get_repulsion_forces()
        
        if self.isSCC == True:
           self.get_scc_forces()

        self.total_forces = [ [0.0 for i in xrange(3)] for j in xrange(self.n_atom) ]

        for iatom in xrange(self.n_atom):
            self.total_forces[iatom][0] = self.electronic_force_x[iatom] + self.repulsion_force_x[iatom]
            self.total_forces[iatom][1] = self.electronic_force_y[iatom] + self.repulsion_force_y[iatom]
            self.total_forces[iatom][2] = self.electronic_force_z[iatom] + self.repulsion_force_z[iatom]

            if self.isSCC == True:

                self.total_forces[iatom][0] = self.total_forces[iatom][0] + self.scc_force_x[iatom]
                self.total_forces[iatom][1] = self.total_forces[iatom][1] + self.scc_force_y[iatom]
                self.total_forces[iatom][2] = self.total_forces[iatom][2] + self.scc_force_z[iatom]

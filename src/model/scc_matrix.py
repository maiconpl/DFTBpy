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

import math
from view.parser import *

import os
_path = os.path.dirname(os.path.realpath(__file__)) + "/" + "../lib/build_scc_matrix.so"
_path2 = os.path.dirname(os.path.realpath(__file__)) + "/" + "../lib/build_gamma_matrix.so"

from ctypes import cdll, byref, c_int, c_double, c_char
lib_scc_matrix = cdll.LoadLibrary(_path)
lib_gamma_matrix = cdll.LoadLibrary(_path2)

class SCCMatrix(object):
      
    def __init__(self, slako_read_list, 
                 atom_type,
                 atoms_type_combinations_list,
                 n_atom,
                 xyz_atom_symbols, 
                 dist_matrix,
                 total_orb, 
                 atom_numb, 
                 H_dftb, 
                 S_dftb,
                 net_charge):

        self.slako_read_list = slako_read_list
        self.atom_type = atom_type
        self.atoms_type_combinations_list = atoms_type_combinations_list
        self.n_atom = n_atom
        self.xyz_atom_symbols = xyz_atom_symbols
        self.dist_matrix = dist_matrix
        self.total_orb = total_orb
        self.atom_numb = atom_numb
        self.H_dftb = H_dftb
        self.S_dftb = S_dftb
        self.net_charge = net_charge

        self.hubbard_s = []
        self.gamma = []
        self.scc_matrix = []
        self.scc_matrix_without_S = []
        self.scc_dftb_matrix = []
        
        # Defining the ctypes

        self._n_atom = c_int(n_atom)
        self._total_orb = c_int(total_orb)
        
        _atom_numb = c_int*(total_orb)
        self._atom_numb = _atom_numb()
        
        _S_dftb = (c_double*(total_orb))*total_orb        
        self._S_dftb = _S_dftb()
        _net_charge = c_double*(n_atom)
        self._net_charge = _net_charge()
        _gamma = (c_double*(n_atom))*n_atom
        self._gamma = _gamma()
        _scc_matrix = (c_double*(total_orb))*total_orb
        self._scc_matrix = _scc_matrix()
        _scc_matrix_without_S = (c_double*(total_orb))*total_orb
        self._scc_matrix_without_S = _scc_matrix_without_S()
        
        _dist_matrix = (c_double*(n_atom))*n_atom
        self._dist_matrix = _dist_matrix()     
        _xyz_atom_symbols = ((c_char*3)*n_atom)
        self._xyz_atom_symbols = _xyz_atom_symbols()
        _hubbard_s = c_double*(n_atom)
        self._hubbard_s = _hubbard_s()
                
        for i in xrange(total_orb):
            self._atom_numb[i] = atom_numb[i]
            for j in xrange(i, total_orb):
                self._S_dftb[i][j] = S_dftb[i][j]
                self._S_dftb[j][i] = S_dftb[i][j]
                self._scc_matrix[i][j] = 0.0
                self._scc_matrix[j][i] = 0.0
                self._scc_matrix_without_S[i][j] = 0.0
                self._scc_matrix_without_S[j][i] = 0.0
                
        for i in xrange(n_atom):
            self._net_charge[i] = net_charge[i]
            self._hubbard_s[i] = 0.0
            self._xyz_atom_symbols[i].value = xyz_atom_symbols[i]
            for j in xrange(i, n_atom):
                self._gamma[i][j] = 0.0
                self._gamma[j][i] = 0.0
                self._dist_matrix[i][j] = dist_matrix[i][j]
                self._dist_matrix[j][i] = dist_matrix[i][j]
                
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

    def build_gamma_matrix_from_C(self):

        for i in xrange(self.n_atom):
            self._hubbard_s[i] = self.hubbard_s[i]

        self.gamma = [ [0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.

        lib_gamma_matrix.build_gamma_matrix(self._n_atom,
                                   byref(self._dist_matrix),
                                   byref(self._xyz_atom_symbols),
                                   byref(self._hubbard_s),
                                   byref(self._gamma))

        self.gamma = [ [0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.

        for i in xrange(self.n_atom):
            for j in xrange(i, self.n_atom):
                self.gamma[i][j] = self._gamma[i][j]
                self.gamma[j][i] = self._gamma[i][j]
                
    def build_gamma_matrix(self):

        self.gamma = [ [0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.
        dist_matrix = [ [0 for i in xrange(self.n_atom)] for j in xrange(self.n_atom)] # defining the dimension of the matrix.
        fac = 0.0

        for in_atom in xrange(self.n_atom):
            for jn_atom in xrange(in_atom, self.n_atom):

                dist_matrix[in_atom][jn_atom] = self.dist_matrix[in_atom][jn_atom]

                ta = 3.2*self.hubbard_s[in_atom]
                tb = 3.2*self.hubbard_s[jn_atom]

                if dist_matrix[in_atom][jn_atom] == 0.0: # the same atom
                   self.gamma[in_atom][jn_atom] = 0.5*( (ta*tb/(ta + tb)) + ( (ta*tb)**2/(ta + tb)**3) )

                if dist_matrix[in_atom][jn_atom] != 0.0 and abs(ta - tb) < 10.0**(-4): # gamma's are very close, e. g. for the same atom type.
                   fac = ( (1.6*dist_matrix[in_atom][jn_atom]*ta*tb)/(ta + tb) )*(1.0 + (ta*tb)/(ta + tb)**2 )
                   self.gamma[in_atom][jn_atom] = 1.0/dist_matrix[in_atom][jn_atom] - (48 + 33*fac + (9.0 + fac)*fac**2)*exp(-fac)/(48*dist_matrix[in_atom][jn_atom])

                if self.xyz_atom_symbols[in_atom] != self.xyz_atom_symbols[jn_atom]: # gamma's are different

                   self.gamma[in_atom][jn_atom] = 1.0/dist_matrix[in_atom][jn_atom] - \
                                                        (math.exp(-(ta*dist_matrix[in_atom][jn_atom]))* \
                                                         ( ((tb**4*ta)/(2*(ta**2 - tb**2)**2)) - \
                                                           ((tb**6 - 3*tb**4*ta**2)/((ta**2 - tb**2)**3*dist_matrix[in_atom][jn_atom])) ))  - \
                                                        (math.exp(-(tb*dist_matrix[in_atom][jn_atom]))* \
                                                         ( ((ta**4*tb)/(2*(tb**2 - ta**2)**2)) - \
                                                           ((ta**6 - 3*ta**4*tb**2)/((tb**2 - ta**2)**3*dist_matrix[in_atom][jn_atom])) ))

                self.gamma[jn_atom][in_atom] = self.gamma[in_atom][jn_atom]

    def build_scc_matrix_from_C(self):
        
        for i in xrange(self.n_atom):
            for j in xrange(self.n_atom):
                self._gamma[i][j] = self.gamma[i][j]
        
        lib_scc_matrix.build_scc_matrix(self._n_atom,
                                   self._total_orb,
                                   byref(self._atom_numb),
                                   byref(self._S_dftb),
                                   byref(self._net_charge),
                                   byref(self._gamma),
                                   byref(self._scc_matrix),
                                   byref(self._scc_matrix_without_S))

        self.scc_matrix = [ [0.0 for i in xrange(len(self.atom_numb))] for j in xrange(len(self.atom_numb))]
        self.scc_matrix_without_S = [ [0.0 for i in xrange(len(self.atom_numb))] for j in xrange(len(self.atom_numb))]

        for i in xrange(self.total_orb):
            for j in xrange(i, self.total_orb):
                self.scc_matrix[i][j] = self._scc_matrix[i][j]
                self.scc_matrix[j][i] = self._scc_matrix[i][j]
                self.scc_matrix_without_S[i][j] = self._scc_matrix_without_S[i][j]
                self.scc_matrix_without_S[j][i] = self._scc_matrix_without_S[i][j]
                        
    def build_scc_matrix(self):

        self.scc_matrix = [ [0 for i in xrange(len(self.atom_numb))] for j in xrange(len(self.atom_numb))] # defining the dimension of the matrix.

        for itotal_orb in xrange(self.total_orb):
            for jtotal_orb in xrange(itotal_orb, self.total_orb):
                for kn_atom in xrange(self.n_atom):
                    self.scc_matrix[itotal_orb][jtotal_orb] = self.scc_matrix[itotal_orb][jtotal_orb] + 0.5*self.S_dftb[itotal_orb][jtotal_orb]*(self.gamma[self.atom_numb[itotal_orb]][kn_atom] + self.gamma[self.atom_numb[jtotal_orb]][kn_atom])*(-self.net_charge[kn_atom])
                self.scc_matrix[jtotal_orb][itotal_orb] = self.scc_matrix[itotal_orb][jtotal_orb]

    def build_scc_dftb_matrix(self):

        self.scc_dftb_matrix = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)] # defining the dimension of the matrix

        for itotal_orb in xrange(self.total_orb):
            for jtotal_orb in xrange(itotal_orb, self.total_orb):
                self.scc_dftb_matrix[itotal_orb][jtotal_orb] = self.H_dftb[itotal_orb][jtotal_orb] +  self.scc_matrix[itotal_orb][jtotal_orb]
                self.scc_dftb_matrix[jtotal_orb][itotal_orb] = self.scc_dftb_matrix[itotal_orb][jtotal_orb]

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
_path = os.path.dirname(os.path.realpath(__file__)) + "/" + "../lib/build_P_matrix.so"

from ctypes import cdll, byref, c_int, c_double
lib_P_matrix = cdll.LoadLibrary(_path)

class DensityMatrix(object):

    def __init__(self, eigen_vectors, eigen_values, total_orb, n_electrons):
  
        self.eigen_vectors = eigen_vectors
        self.eigen_values = eigen_values
        self.total_orb = total_orb
        self.n_electrons = n_electrons
        self.P_dftb = []
        self.P_dftb_e = []
        
        # Defining the ctypes

        self._n_electrons = c_int(int(n_electrons))
        self._total_orb = c_int(total_orb)

        _eigen_vectors = (c_double*(total_orb))*total_orb
        self._eigen_vectors = _eigen_vectors()

        _eigen_values = c_double*(total_orb)
        self._eigen_values = _eigen_values()
        
        _P_dftb = (c_double*(total_orb))*total_orb  
        self._P_dftb = _P_dftb()

        _P_dftb_e = (c_double*(total_orb))*total_orb         
        self._P_dftb_e = _P_dftb_e()
               
        for i in xrange(total_orb):
            self._eigen_values[i] = eigen_values[i]
            for j in xrange(total_orb):
                self._P_dftb[i][j] = 0.0
                self._P_dftb_e[i][j] = 0.0
                self._eigen_vectors[i][j] = eigen_vectors[i][j]
    
    def build_P_matrix_from_C(self):

        lib_P_matrix.build_p_matrix(self._n_electrons,
                                 self._total_orb,
                                 byref(self._eigen_vectors),
                                 byref(self._eigen_values),
                                 byref(self._P_dftb),
                                 byref(self._P_dftb_e))

        self.P_dftb = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]
        self.P_dftb_e = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        for i in xrange(self.total_orb):
            for j in xrange(i, self.total_orb):
                self.P_dftb[i][j] = self._P_dftb[i][j]
                self.P_dftb[j][i] = self.P_dftb[i][j]
                self.P_dftb_e[i][j] = self._P_dftb_e[i][j]
                self.P_dftb_e[j][i] = self.P_dftb_e[i][j]
                
    def build_P_matrix(self):

        self.P_dftb = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)]

        # Getting the right (upper) triangle of P matrix
        for i in xrange(self.total_orb):
            for j in xrange(i, self.total_orb):
                self.P_dftb[i][j] = 0.0
                for k in xrange(int(self.n_electrons/2)): # closed shell
                    self.P_dftb[i][j] = self.P_dftb[i][j] + 2.0*self.eigen_vectors[i][k]*self.eigen_vectors[j][k]
                self.P_dftb[j][i] = self.P_dftb[i][j] # Getting the left (lower) triangle of P matrix

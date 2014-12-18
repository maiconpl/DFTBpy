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

from numpy import *
from scipy import linalg

# External export, from fortran or other languages
#from lib import get_diag

class Diagonalization(object):
    """ This class is responsable to the diagonalization of the H and S matrix"""

    def __init__(self, H_dftb, S_dftb, total_orb):

        self.H_dftb = H_dftb
        self.S_dftb = S_dftb
        self.total_orb = total_orb

        self.eigen_vectors = []
        self.eigen_values = []

    def dftb_matrix_diagonalization(self):

        self.H_dftb = array(self.H_dftb)
        self.S_dftb = array(self.S_dftb)
        self.eigen_vectors = array(self.eigen_vectors)

#        eigen_vectors = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)] # defining the dimension of the matrix.
#        eigen_values  = [0 for i in xrange(self.total_orb)] # defining the dimension of the matrix.

#        self.eigen_vectors, self.eigen_values = get_diag.get_dsygvd(self.H_dftb, self.S_dftb, eigen_vectors, eigen_values, self.total_orb)
        v, w  = linalg.eigh(self.H_dftb, self.S_dftb)
        self.eigen_values = v.real
        self.eigen_vectors = w.real

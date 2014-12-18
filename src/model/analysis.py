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

import numpy as np

class Analysis(object):

    def __init__(self, atom_type, 
                 xyz_atom_symbols, 
                 max_angular_momentum, 
                 total_orb, 
                 n_electrons_atom, 
                 S_dftb, 
                 P_dftb):
  
        self.atom_type = atom_type
        self.xyz_atom_symbols = xyz_atom_symbols
        self.max_angular_momentum = max_angular_momentum
        self.total_orb = total_orb
        self.n_electrons_atom = n_electrons_atom
        self.S_dftb = S_dftb
        self.P_dftb = P_dftb

        self.PS = []
        self.N_Mulliken = []
        self.net_charge = []

    def mulliken_analysis(self):

        self.PS = np.dot(self.P_dftb,self.S_dftb) # that's the multiplication of the P and S matrix
        self.N_Mulliken = [0 for i in xrange(len(self.xyz_atom_symbols))]
        self.net_charge = [0 for i in xrange(len(self.xyz_atom_symbols))]
 
        n_orb = 0 # number of orbital per atom
        total_n_orb = 0 # total orbital per atom
        tmp1 = 0 # temporary variable
        tmp2 = 0.0 # temporary variable

        for ixyz_atom_symbols in xrange(len(self.xyz_atom_symbols)):
            for jatom_type in xrange(len(self.atom_type)):
                if self.xyz_atom_symbols[ixyz_atom_symbols] == self.atom_type[jatom_type]:
                   if self.max_angular_momentum[jatom_type] == "s":
                      n_orb = 1
                   if self.max_angular_momentum[jatom_type] == "p":
                      n_orb = 4
                   if self.max_angular_momentum[jatom_type] == "d":
                      n_orb = 9
                   total_n_orb = total_n_orb + n_orb
                   for itotal_n_orb in xrange(tmp1, total_n_orb):
                       self.N_Mulliken[ixyz_atom_symbols] = tmp2 + self.PS[itotal_n_orb][itotal_n_orb]
                       self.net_charge[ixyz_atom_symbols] = self.n_electrons_atom[ixyz_atom_symbols] - self.N_Mulliken[ixyz_atom_symbols]
                       tmp2 = self.N_Mulliken[ixyz_atom_symbols]
                   tmp2 = 0
                   tmp1 = total_n_orb

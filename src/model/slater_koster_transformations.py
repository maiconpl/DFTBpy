"""
Created on Oct 25, 2014 by MPL.

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

from model.fit import Fit
from math import sqrt

class SlaterKosterTransformations(object):

    def __init__(self, dist_list,
                 dist_step_list,
                 H_list,
                 S_list,
                 max_all_momentum):

        self.dist_list = dist_list
        self.dist_step_list = dist_step_list
        self.H_list = H_list
        self.S_list = S_list
        self.max_all_momentum = max_all_momentum

        self.fit_slako = Fit()
    
    def get_skf(self, atom_type_numb_A,
                      atom_type_numb_B,
                      integral_symmetry, 
                      dist_atoms,
                      dist_x,
                      dist_y,
                      dist_z):
        
        H_dftb = 0.0
        S_dftb = 0.0
        
        if self.max_all_momentum >= 0:
                      
           # ss 1
           if integral_symmetry == "ss":
              #
              _H_9 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][9],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][9],
                                                               dist_atoms)

              _S_9 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][9],
                                                            self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                            self.S_list[atom_type_numb_A][atom_type_numb_B][9],
                                                            dist_atoms)
              #
              H_dftb = _H_9
              S_dftb = _S_9
              
              return H_dftb, S_dftb

        if self.max_all_momentum >= 1:

           # spx 2
           if integral_symmetry == "spx":
              #
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)
              #
              H_dftb = -(dist_x/dist_atoms)*(_H_8)
              S_dftb = -(dist_x/dist_atoms)*(_S_8)
              
              return H_dftb, S_dftb

           # pxs 3
           if integral_symmetry == "pxs":
              # Inverting the orbitals
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)
              #
              H_dftb = (dist_x/dist_atoms)*(_H_8)
              S_dftb = (dist_x/dist_atoms)*(_S_8)

              return H_dftb, S_dftb

           # spy 4
           if integral_symmetry == "spy":
              #
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)
              #
              H_dftb = -(dist_y/dist_atoms)*(_H_8)
              S_dftb = -(dist_y/dist_atoms)*(_S_8)

              return H_dftb, S_dftb
   
           # pys 5
           if integral_symmetry == "pys":
              # Inverting the orbitals
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)
              #
              H_dftb = (dist_y/dist_atoms)*(_H_8)
              S_dftb = (dist_y/dist_atoms)*(_S_8)

              return H_dftb, S_dftb

           # spz 6
           if integral_symmetry == "spz":
              #
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][8],
                                                               dist_atoms)
              #
              H_dftb = -(dist_z/dist_atoms)*(_H_8)
              S_dftb = -(dist_z/dist_atoms)*(_S_8)

              return H_dftb, S_dftb

           # pzs 7
           if integral_symmetry == "pzs":
              # Inverting the orbitals
              _H_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)

              _S_8 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][8],
                                                               dist_atoms)
              #
              H_dftb = (dist_z/dist_atoms)*(_H_8)
              S_dftb = (dist_z/dist_atoms)*(_S_8)

              return H_dftb, S_dftb

           # pxpx 8
           if integral_symmetry == "pxpx":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	      #
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = ((dist_x/dist_atoms)**2)*(_H_5) +\
                        (1 - (dist_x/dist_atoms)**2)*(_H_6)
              S_dftb = ((dist_x/dist_atoms)**2)*(_S_5) +\
                        (1 - (dist_x/dist_atoms)**2)*(_S_6)

              return H_dftb, S_dftb

           # pxpy 9
           if integral_symmetry == "pxpy":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	      #
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_x/dist_atoms)*(dist_y/dist_atoms)*(_H_5) -\
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)*(_H_6)
              S_dftb = (dist_x/dist_atoms)*(dist_y/dist_atoms)*(_S_5) -\
                                           (dist_x/dist_atoms)*(dist_y/dist_atoms)*(_S_6)

              return H_dftb, S_dftb

           # pypx 10
           if integral_symmetry == "pypx":
              # Inverting the orbitals
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
	      #
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_y/dist_atoms)*(dist_x/dist_atoms)*(_H_5) -\
                       (dist_y/dist_atoms)*(dist_x/dist_atoms)*(_H_6)
              S_dftb = (dist_y/dist_atoms)*(dist_x/dist_atoms)*(_S_5) -\
                       (dist_y/dist_atoms)*(dist_x/dist_atoms)*(_S_6)

              return H_dftb, S_dftb

           # pxpz 11
           if integral_symmetry == "pxpz":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_x/dist_atoms)*(dist_z/dist_atoms)*(_H_5) -\
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)*(_H_6)
              S_dftb = (dist_x/dist_atoms)*(dist_z/dist_atoms)*(_S_5) -\
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)*(_S_6)

              return H_dftb, S_dftb

           # pzpx 12
           if integral_symmetry == "pzpx":
              # Inverting the orbitals
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_z/dist_atoms)*(dist_x/dist_atoms)*(_H_5) -\
                       (dist_z/dist_atoms)*(dist_x/dist_atoms)*(_H_6)
              S_dftb = (dist_z/dist_atoms)*(dist_x/dist_atoms)*(_S_5) -\
                       (dist_z/dist_atoms)*(dist_x/dist_atoms)*(_S_6)

              return H_dftb, S_dftb

           # pypy 13
           if integral_symmetry == "pypy":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = ((dist_y/dist_atoms)**2)*(_H_5) +\
                        (1 - (dist_y/dist_atoms)**2)*(_H_6)
              S_dftb = ((dist_y/dist_atoms)**2)*(_S_5) +\
                        (1 - (dist_y/dist_atoms)**2)*(_S_6)

              return H_dftb, S_dftb

           # pypz 14
           if integral_symmetry == "pypz":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_y/dist_atoms)*(dist_z/dist_atoms)*(_H_5) -\
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)*(_H_6)
              S_dftb = (dist_y/dist_atoms)*(dist_z/dist_atoms)*(_S_5) -\
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)*(_S_6)

              return H_dftb, S_dftb

           # pzpy 15
           if integral_symmetry == "pzpy":
              # Inverting the orbitals
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][6],
                                                               dist_atoms)
              #
              H_dftb = (dist_z/dist_atoms)*(dist_y/dist_atoms)*(_H_5) -\
                       (dist_z/dist_atoms)*(dist_y/dist_atoms)*(_H_6)
              S_dftb = (dist_z/dist_atoms)*(dist_y/dist_atoms)*(_S_5) -\
                       (dist_z/dist_atoms)*(dist_y/dist_atoms)*(_S_6)

              return H_dftb, S_dftb
   
           # pzpz 16
           if integral_symmetry == "pzpz":
              #
              _H_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)

              _S_5 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][5],
                                                               dist_atoms)
               
              _H_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
	       
              _S_6 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][6],
                                                               dist_atoms)
              #
              H_dftb = ((dist_z/dist_atoms)**2)*(_H_5) +\
                        (1 - (dist_z/dist_atoms)**2)*(_H_6)
              S_dftb = ((dist_z/dist_atoms)**2)*(_S_5) +\
                        (1 - (dist_z/dist_atoms)**2)*(_S_6)

              return H_dftb, S_dftb

        if self.max_all_momentum == 2:
           # sdxy 17
           if integral_symmetry == "sdxy":
              #
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_y/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_y/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # dxys 18
           if integral_symmetry == "dxys":
              # Inverting the orbitals
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_y/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_y/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # sdxz 19
           if integral_symmetry == "sdxz":
              #
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # dxzs 20
           if integral_symmetry == "dxzs":
              # Inverting the orbitals
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # sdyz 21
           if integral_symmetry == "sdyz":
              #
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_y/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_y/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # dyzs 22
           if integral_symmetry == "dyzs":
              # Inverting the orbitals
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_y/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_H_7)
              S_dftb = sqrt(3)*(dist_y/dist_atoms)* \
                               (dist_z/dist_atoms)* \
                               (_S_7)

              return H_dftb, S_dftb

           # sdx2-y2 23
           if integral_symmetry == "sdx2-y2":
              #
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)* \
                       ( (dist_x/dist_atoms)**2 - \
                       (dist_y/dist_atoms)**2 )* \
                       (_H_7)
              S_dftb = 0.5*sqrt(3)* \
                       ( (dist_x/dist_atoms)**2 - \
                       (dist_y/dist_atoms)**2 )* \
                       (_S_7)

              return H_dftb, S_dftb

           # dx2-y2s 24
           if integral_symmetry == "dx2-y2s":
              # Inverting the orbitals
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)* \
                       ( (dist_x/dist_atoms)**2 - \
                       (dist_y/dist_atoms)**2 )* \
                       (_H_7)
              S_dftb = 0.5*sqrt(3)* \
                       ( (dist_x/dist_atoms)**2 - \
                       (dist_y/dist_atoms)**2 )* \
                       (_S_7)

              return H_dftb, S_dftb

           # sdz2 25
           if integral_symmetry == "sdz2":
              #
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][7],
                                                               dist_atoms)
              #
              H_dftb = ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + \
                       (dist_y/dist_atoms)**2 ) )* \
                       (_H_7)
              S_dftb = ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + \
                       (dist_y/dist_atoms)**2 ) )* \
                       (_S_7)

              return H_dftb, S_dftb

           # dz2s 26
           if integral_symmetry == "dz2s":
              # Inverting the orbitals
              _H_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
	       
              _S_7 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][7],
                                                               dist_atoms)
              #
              H_dftb = ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + \
                       (dist_y/dist_atoms)**2 ) )* \
                       (_H_7)
              S_dftb = ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + \
                       (dist_y/dist_atoms)**2 ) )* \
                       (_S_7)

              return H_dftb, S_dftb

           # pxdxy 27
           if integral_symmetry == "pxdxy":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*( (dist_x/dist_atoms)**2 * \
                       (dist_y/dist_atoms) )* \
                       (_H_3) - \
                       (dist_y/dist_atoms)* \
                       ( 1 - 2*( (dist_x/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = -sqrt(3)*( (dist_x/dist_atoms)**2 * \
                       (dist_y/dist_atoms) )* \
                       (_S_3) - \
                       (dist_y/dist_atoms)* \
                       ( 1 - 2*( (dist_x/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxypx 28
           if integral_symmetry == "dxypx":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*( (dist_x/dist_atoms)**2 * \
                       (dist_y/dist_atoms) )* \
                       (_H_3) + \
                       (dist_y/dist_atoms)* \
                       ( 1 - 2*( (dist_x/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = sqrt(3)*( (dist_x/dist_atoms)**2 * \
                       (dist_y/dist_atoms) )* \
                       (_S_3) + \
                       (dist_y/dist_atoms)* \
                       ( 1 - 2*( (dist_x/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pydxy 29
           if integral_symmetry == "pydxy":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*((dist_y/dist_atoms)**2)*(dist_x/dist_atoms)* \
                       (_H_3) - \
                       (dist_x/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = -sqrt(3)*((dist_y/dist_atoms)**2)*(dist_x/dist_atoms)* \
                       (_S_3) - \
                       (dist_x/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxypy 30
           if integral_symmetry == "dxypy":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*((dist_y/dist_atoms)**2)*(dist_x/dist_atoms)* \
                       (_H_3) + \
                       (dist_x/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = sqrt(3)*((dist_y/dist_atoms)**2)*(dist_x/dist_atoms)* \
                       (_S_3) + \
                       (dist_x/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pzdxy 31
           if integral_symmetry == "pzdxy":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*( (dist_x/dist_atoms)* \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms) )* \
                       (_H_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       (dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = -sqrt(3)*( (dist_x/dist_atoms)* \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms) )* \
                       (_S_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       (dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxypz 32
           if integral_symmetry == "dxypz":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*( (dist_x/dist_atoms)* \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms) )* \
                       (_H_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       (dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = sqrt(3)*( (dist_x/dist_atoms)* \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms) )* \
                       (_S_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       (dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pxdyz 33
           if integral_symmetry == "pxdyz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = -sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dyzpx 34
           if integral_symmetry == "dyzpx":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pydyz 35
           if integral_symmetry == "pydyz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_3) - \
                       (dist_z/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = -sqrt(3)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_3) - \
                       (dist_z/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dyzpy 36
           if integral_symmetry == "dyzpy":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_3) + \
                       (dist_z/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = sqrt(3)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_3) + \
                       (dist_z/dist_atoms)* \
                       (1 - 2*((dist_y/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pzdyz 37
           if integral_symmetry == "pzdyz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*((dist_z/dist_atoms)**2)*(dist_y/dist_atoms)* \
                       (_H_3) - \
                       (dist_y/dist_atoms)* \
                       (1 - 2*((dist_z/dist_atoms)**2) )* \
                       (_H_4) 
              S_dftb = -sqrt(3)*((dist_z/dist_atoms)**2)*(dist_y/dist_atoms)* \
                       (_S_3) - \
                       (dist_y/dist_atoms)* \
                       (1 - 2*((dist_z/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dyzpz 38
           if integral_symmetry == "dyzpz":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*((dist_z/dist_atoms)**2)*(dist_y/dist_atoms)* \
                       (_H_3) + \
                       (dist_y/dist_atoms)* \
                       (1 - 2*((dist_z/dist_atoms)**2) )* \
                       (_H_4) 
              S_dftb = sqrt(3)*((dist_z/dist_atoms)**2)*(dist_y/dist_atoms)* \
                       (_S_3) + \
                       (dist_y/dist_atoms)* \
                       (1 - 2*((dist_z/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pxdxz 39
           if integral_symmetry == "pxdxz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*((dist_x/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_3) - \
                       (dist_z/dist_atoms)* \
                       ( 1 - 2*((dist_x/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = -sqrt(3)*((dist_x/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_3) - \
                       (dist_z/dist_atoms)* \
                       ( 1 - 2*((dist_x/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxzpx 40
           if integral_symmetry == "dxzpx":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*((dist_x/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_3) + \
                       (dist_z/dist_atoms)* \
                       ( 1 - 2*((dist_x/dist_atoms)**2) )* \
                       (_H_4)
              S_dftb = sqrt(3)*((dist_x/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_3) + \
                       (dist_z/dist_atoms)* \
                       ( 1 - 2*((dist_x/dist_atoms)**2) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pydxz 41
           if integral_symmetry == "pydxz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = -sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_3) + \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxzpy 42
           if integral_symmetry == "dxzpy":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_4)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_3) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pzdxz 43
           if integral_symmetry == "pzdxz":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_3) - \
                       (dist_x/dist_atoms)* \
                       ( 1 - 2*(dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = -sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_3) - \
                       (dist_x/dist_atoms)* \
                       ( 1 - 2*(dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxzpz 44
           if integral_symmetry == "dxzpz":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_3) + \
                       (dist_x/dist_atoms)* \
                       ( 1 - 2*(dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_3) + \
                       (dist_x/dist_atoms)* \
                       ( 1 - 2*(dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pxdx2-y2 45
           if integral_symmetry == "pxdx2-y2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -0.5*sqrt(3)*(dist_x/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) - \
                       (dist_x/dist_atoms)* \
                       (1 - ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_H_4)
              S_dftb = -0.5*sqrt(3)*(dist_x/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) - \
                       (dist_x/dist_atoms)* \
                       (1 - ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dx2-y2px 46
           if integral_symmetry == "dx2-y2px":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)*(dist_x/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) + \
                       (dist_x/dist_atoms)* \
                       (1 - ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_H_4)
              S_dftb = 0.5*sqrt(3)*(dist_x/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) + \
                       (dist_x/dist_atoms)* \
                       (1 - ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pydx2-y2 47
           if integral_symmetry == "pydx2-y2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -0.5*sqrt(3)*(dist_y/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) + \
                       (dist_y/dist_atoms)* \
                       (1 + ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_H_4)
              S_dftb = -0.5*sqrt(3)*(dist_y/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) + \
                       (dist_y/dist_atoms)* \
                       (1 + ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dx2-y2py 48
           if integral_symmetry == "dx2-y2py":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)*(dist_y/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) - \
                       (dist_y/dist_atoms)* \
                       (1 + ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_H_4)
              S_dftb = 0.5*sqrt(3)*(dist_y/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) - \
                       (dist_y/dist_atoms)* \
                       (1 + ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pzdx2-y2 49
           if integral_symmetry == "pzdx2-y2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -0.5*sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) + \
                       (dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) * \
                       (_H_4)
              S_dftb = -0.5*sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) + \
                       (dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dx2-y2pz 50
           if integral_symmetry == "dx2-y2pz":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_H_3) - \
                       (dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 ) * \
                       (_H_4)
              S_dftb = 0.5*sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_3) - \
                       (dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 - (dist_y/dist_atoms)**2 )* \
                       (_S_4)

              return H_dftb, S_dftb

           # pxdz2 51
           if integral_symmetry == "pxdz2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -(dist_x/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) + \
                       sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = -(dist_x/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) + \
                       sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dz2px 52
           if integral_symmetry == "dz2px":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = (dist_x/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) - \
                       sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = (dist_x/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) - \
                       sqrt(3)*(dist_x/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pydz2 53
           if integral_symmetry == "pydz2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -(dist_y/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) + \
                       sqrt(3)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = -(dist_y/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) + \
                       sqrt(3)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # dz2py 54
           if integral_symmetry == "dz2py":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = (dist_y/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) - \
                       sqrt(3)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_4)
              S_dftb = (dist_y/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) - \
                       sqrt(3)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_4)

              return H_dftb, S_dftb

           # pzdz2 55
           if integral_symmetry == "pzdz2":
              #
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][4],
                                                               dist_atoms)
              #
              H_dftb = -(dist_z/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) - \
                       sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 )* \
                       (_H_4)
              S_dftb = -(dist_z/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) - \
                       sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dz2pz 56
           if integral_symmetry == "dz2pz":
              # Inverting the orbitals
              _H_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
	       
              _S_3 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][3],
                                                               dist_atoms)
               
              _H_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
	       
              _S_4 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][4],
                                                               dist_atoms)
              #
              H_dftb = (dist_z/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_H_3) + \
                       sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 )* \
                       (_H_4)
              S_dftb = (dist_z/dist_atoms)* \
                       ( (dist_z/dist_atoms)**2 - \
                       0.5*( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) )* \
                       (_S_3) + \
                       sqrt(3)*(dist_z/dist_atoms)* \
                       ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 )* \
                       (_S_4)

              return H_dftb, S_dftb

           # dxydxy 57
           if integral_symmetry == "dxydxy":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2)* \
                       (_H_0) + \
                       ( ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) - \
                       4*((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       ( ((dist_z/dist_atoms)**2) + \
                       ((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 3*((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2)* \
                       (_S_0) + \
                       ( ( (dist_x/dist_atoms)**2 + (dist_y/dist_atoms)**2 ) - \
                       4*((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       ( ((dist_z/dist_atoms)**2) + \
                       ((dist_x/dist_atoms)**2)*((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dxydyz 58
           if integral_symmetry == "dxydyz":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*(dist_x/dist_atoms)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_y/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_y/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*(dist_x/dist_atoms)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_y/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_y/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dyzdxy 59
           if integral_symmetry == "dyzdxy":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 3*(dist_x/dist_atoms)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_y/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_y/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*(dist_x/dist_atoms)*((dist_y/dist_atoms)**2)*(dist_z/dist_atoms)* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_y/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_y/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dxydxz 60
           if integral_symmetry == "dxydxz":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*((dist_x/dist_atoms)**2)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_0) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_x/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*((dist_x/dist_atoms)**2)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_0) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_x/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dxzdxy 61
           if integral_symmetry == "dxzdxy":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 3*((dist_x/dist_atoms)**2)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_H_0) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_x/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*((dist_x/dist_atoms)**2)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       (_S_0) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 4*(dist_x/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dyzdyz 62
           if integral_symmetry == "dyzdyz":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2)* \
                       (_H_0) + \
                       ( ( (dist_y/dist_atoms)**2 + (dist_z/dist_atoms)**2 ) - \
                       4*((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_H_1) + \
                       ( ((dist_x/dist_atoms)**2) + \
                       ((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 3*((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2)* \
                       (_S_0) + \
                       ( ( (dist_y/dist_atoms)**2 + (dist_z/dist_atoms)**2 ) - \
                       4*((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_S_1) + \
                       ( ((dist_x/dist_atoms)**2) + \
                       ((dist_y/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dyzdxz 63
           if integral_symmetry == "dyzdxz":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*(dist_x/dist_atoms)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 - 4*(dist_z/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*(dist_x/dist_atoms)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 - 4*(dist_z/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dxzdyz 64
           if integral_symmetry == "dxzdyz":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 3*(dist_x/dist_atoms)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 - 4*(dist_z/dist_atoms)**2 )* \
                       (_H_1) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - 1 ) * \
                       (_H_2)
              S_dftb = 3*(dist_x/dist_atoms)*(dist_y/dist_atoms)*((dist_z/dist_atoms)**2)* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 - 4*(dist_z/dist_atoms)**2 )* \
                       (_S_1) + \
                       (dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - 1 ) * \
                       (_S_2)

              return H_dftb, S_dftb

           # dxzdxz 65
           if integral_symmetry == "dxzdxz":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 3*((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2)* \
                       (_H_0) + \
                       ( ( (dist_x/dist_atoms)**2 + (dist_z/dist_atoms)**2 ) - \
                       4*((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_H_1) + \
                       ( ((dist_y/dist_atoms)**2) + \
                       ((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 3*((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2)* \
                       (_S_0) + \
                       ( ( (dist_x/dist_atoms)**2 + (dist_z/dist_atoms)**2 ) - \
                       4*((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_S_1) + \
                       ( ((dist_y/dist_atoms)**2) + \
                       ((dist_x/dist_atoms)**2)*((dist_z/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dxydx2-y2 66
           if integral_symmetry == "dxydx2-y2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       0.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       0.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dx2-y2dxy 67
           if integral_symmetry == "dx2-y2dxy":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       0.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) - \
                       2*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       0.5*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dyzdx2-y2 68
           if integral_symmetry == "dyzdx2-y2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) - \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) - \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dx2-y2dyz 69
           if integral_symmetry == "dx2-y2dyz":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) - \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) - \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_1) + \
                       (dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 + 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dxzdx2-y2 70
           if integral_symmetry == "dxzdx2-y2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_1) - \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_1) - \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dx2-y2dxz 71
           if integral_symmetry == "dx2-y2dxz":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 1.5*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_1) - \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_H_2)
              S_dftb = 1.5*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_0) + \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 2*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_1) - \
                       (dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( 1 - 0.5*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) ) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dxydz2 72
           if integral_symmetry == "dxydz2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) - \
                       2*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ((dist_z/dist_atoms)**2)*(_H_1) + \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 + ((dist_z/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) - \
                       2*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ((dist_z/dist_atoms)**2)*(_S_1) + \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 + ((dist_z/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dz2dxy 73
           if integral_symmetry == "dz2dxy":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) - \
                       2*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ((dist_z/dist_atoms)**2)*(_H_1) + \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 + ((dist_z/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) - \
                       2*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ((dist_z/dist_atoms)**2)*(_S_1) + \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_y/dist_atoms)* \
                       ( 1 + ((dist_z/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dyzdz2 74
           if integral_symmetry == "dyzdz2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) + \
                       sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_H_1) - \
                       0.5*sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) + \
                       sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_S_1) - \
                       0.5*sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dz2dyz 75
           if integral_symmetry == "dz2dyz":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) + \
                       sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_H_1) - \
                       0.5*sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) + \
                       sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_S_1) - \
                       0.5*sqrt(3)*(dist_y/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dxzdz2 76
           if integral_symmetry == "dxzdz2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) + \
                       sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_H_1) - \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) + \
                       sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_S_1) - \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dz2dxz 77
           if integral_symmetry == "dz2dxz":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) + \
                       sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_H_1) - \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_z/dist_atoms)**2) - \
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) + \
                       sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ((dist_z/dist_atoms)**2) )* \
                       (_S_1) - \
                       0.5*sqrt(3)*(dist_x/dist_atoms)*(dist_z/dist_atoms)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dx2-y2dx2-y2 78
           if integral_symmetry == "dx2-y2dx2-y2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 0.75*( ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2)* \
                       (_H_0) + \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2 )* \
                       (_H_1) + \
                       ( ((dist_z/dist_atoms)**2) + \
                       0.25*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2 )* \
                       (_H_2)
              S_dftb = 0.75*( ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2)* \
                       (_S_0) + \
                       ( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) - \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2 )* \
                       (_S_1) + \
                       ( ((dist_z/dist_atoms)**2) + \
                       0.25*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )**2 )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dx2-y2dz2 79
           if integral_symmetry == "dx2-y2dz2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) - \
                       sqrt(3)*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       0.25*sqrt(3)*( 1 + ((dist_z/dist_atoms)**2) )* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 0.5*sqrt(3)*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) - \
                       sqrt(3)*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       0.25*sqrt(3)*( 1 + ((dist_z/dist_atoms)**2) )* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dz2dx2-y2 80
           if integral_symmetry == "dz2dx2-y2":
              # Inverting the orbitals
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.H_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               self.dist_step_list[atom_type_numb_B][atom_type_numb_A],
                                                               self.S_list[atom_type_numb_B][atom_type_numb_A][2],
                                                               dist_atoms)
              #
              H_dftb = 0.5*sqrt(3)*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_H_0) - \
                       sqrt(3)*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       0.25*sqrt(3)*( 1 + ((dist_z/dist_atoms)**2) )* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_H_2)
              S_dftb = 0.5*sqrt(3)*( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )* \
                       (_S_0) - \
                       sqrt(3)*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       0.25*sqrt(3)*( 1 + ((dist_z/dist_atoms)**2) )* \
                       ( ((dist_x/dist_atoms)**2) - ((dist_y/dist_atoms)**2) )* \
                       (_S_2)

              return H_dftb, S_dftb

           # dz2dz2 81
           if integral_symmetry == "dz2dz2":
              #
              _H_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
	       
              _S_0 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][0],
                                                               dist_atoms)
               
              _H_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)
	       
              _S_1 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][1],
                                                               dist_atoms)

              _H_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.H_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
	       
              _S_2 = self.fit_slako.slako_integrals_fit_from_C(self.dist_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               self.dist_step_list[atom_type_numb_A][atom_type_numb_B],
                                                               self.S_list[atom_type_numb_A][atom_type_numb_B][2],
                                                               dist_atoms)
              #
              H_dftb = ( ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )**2)* \
                       (_H_0) + \
                       3*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_H_1) + \
                       0.75*( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )**2)* \
                       (_H_2)
              S_dftb = ( ( ((dist_z/dist_atoms)**2) -
                       0.5*( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) ) )**2)* \
                       (_S_0) + \
                       3*((dist_z/dist_atoms)**2)* \
                       ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )* \
                       (_S_1) + \
                       0.75*( ( ((dist_x/dist_atoms)**2) + ((dist_y/dist_atoms)**2) )**2)* \
                       (_S_2)

              return H_dftb, S_dftb

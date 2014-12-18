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

from view.parser import *
import sys

class Energy(object):

    def __init__(self, slako_read_list,
                 atom_type,
                 atoms_type_combinations_list,
                 eigen_values, 
                 n_electrons, 
                 n_atom,
                 xyz_atom_symbols, 
                 dist_matrix):

        self.slako_read_list = slako_read_list
        self.atom_type = atom_type
        self.atoms_type_combinations_list = atoms_type_combinations_list
        self.eigen_values = eigen_values
        self.n_electrons = n_electrons

        self.n_atom = n_atom
        self.xyz_atom_symbols = xyz_atom_symbols
        self.dist_matrix = dist_matrix

        self.dftb_total_energy = 0.0
        self.band_energy = 0.0
        self.repulsion_energy = 0.0
        self.rep = 0.0
    
    def get_band_energy(self):

        n_ocupied_orbitals = int(self.n_electrons/2)
        n_unocupied_orbitals = len(self.eigen_values) - n_ocupied_orbitals 
        ocupied_orbitals_energies = []

        if int(self.n_electrons/2) > len(self.eigen_values):

           print ""
           print "ERROR!!!"
           print "The number of (valence electrons)/2 " + "= " + str(self.n_electrons/2)
           print "is greater than the total number of valence orbitals " + "'" + str(len(self.eigen_values)) + "'" + "."
           print "Please, check the number of valence electrons"
           print "and the number of valence atomic orbitals present in the SLAKO."
           print "Also, check in the input if you set the correct number of total charge of the system."
           print ""
           sys.exit()

        for iocupied in xrange(n_ocupied_orbitals):
            ocupied_orbitals_energies.append(self.eigen_values[iocupied])

        self.band_energy = 2*sum(ocupied_orbitals_energies)

        unocupied_orbitals_energies = []
        for iunocupied in xrange(n_ocupied_orbitals,  len(self.eigen_values)):
            unocupied_orbitals_energies.append(self.eigen_values[iunocupied])

        return self.band_energy

    def get_repulsion_energy(self):

        n_atom_type = len(self.atom_type)
        self.rep = 0.0
        for iatom in xrange(self.n_atom):
            for jatom in xrange(iatom + 1, self.n_atom):

                poly_rep = 0.0

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
                    poly_rep = slako_read.poly_coef[0]*(diff**2) + slako_read.poly_coef[1]*(diff**3) + \
                               slako_read.poly_coef[2]*(diff**4) + slako_read.poly_coef[3]*(diff**5) + \
                               slako_read.poly_coef[4]*(diff**6) + slako_read.poly_coef[5]*(diff**7) + \
                               slako_read.poly_coef[6]*(diff**8)
                else:
                    poly_rep = 0.0

                self.repulsion_energy = self.repulsion_energy + poly_rep

        self.repulsion_energy = self.repulsion_energy

        return self.repulsion_energy

    def get_dftb_total_energy(self):
 
        self.dftb_total_energy = self.band_energy + self.repulsion_energy
 
        return self.dftb_total_energy
             
# On for debugging purposes
#
#        print "Number of ocupied orbitals:", n_ocupied_orbitals
#        print "Number of unocupied orbitals:", n_unocupied_orbitals
#        print "Ocupied orbitals:", ocupied_orbitals_energies[0:]
#        print "Unocupied orbitals:", unocupied_orbitals_energies[0:]
#        print "DFTB electronic energy", self.band_energy

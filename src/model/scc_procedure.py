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

from diag import Diagonalization
from density import DensityMatrix
from analysis import Analysis
from scc_matrix import SCCMatrix
import math
import numpy as np
import sys

class SCCProcedure(object):
      
    def __init__(self, slako_read_list,
                 atoms_type_combinations_list, 
                 scc_step,
                 scc_tolerance,
                 atom_type,
                 xyz_atom_symbols, 
                 dist_matrix,
                 max_angular_momentum, 
                 n_electrons_atom,
                 scc_dftb_matrix,
                 H_dftb,
                 S_dftb,
                 gamma,
                 total_orb, 
                 atom_numb, 
                 n_atom,
                 n_electrons, 
                 band_energy, 
                 repulsion_energy,
                 scc_mixing_value,
                 debug):

        self.slako_read_list = slako_read_list
        self.atoms_type_combinations_list = atoms_type_combinations_list
        self.scc_step = scc_step
        self.scc_tolerance = scc_tolerance
        self.atom_type = atom_type
        self.xyz_atom_symbols = xyz_atom_symbols
        self.dist_matrix = dist_matrix
        self.max_angular_momentum = max_angular_momentum
        self.n_electrons_atom = n_electrons_atom
        self.scc_dftb_matrix = scc_dftb_matrix
        self.H_dftb = H_dftb
        self.S_dftb = S_dftb
        self.gamma = gamma
        self.total_orb = total_orb
        self.atom_numb = atom_numb
        self.n_atom = n_atom
        self.n_electrons = n_electrons
        self.band_energy = band_energy
        self.repulsion_energy = repulsion_energy
        self.scc_mixing_value = scc_mixing_value
        self.debug = debug

        self.scc_matrix_without_S_final = [ [0.0 for i in xrange(total_orb)] for j in xrange(total_orb)]
        self.P_dftb_final = [ [0.0 for i in xrange(total_orb)] for j in xrange(total_orb)]        
        self.P_dftb_e_final = [ [0.0 for i in xrange(total_orb)] for j in xrange(total_orb)]

        self.H0_energy = 0.0
        self.scc_energy = 0.0
        self.net_charge = 0.0
        
        self.isMixing = True

    def scc_procedure(self):

        eigen_values_final  = [0 for i in xrange(self.total_orb)] # defining the dimension of the matrix.
        eigen_vectors_final = [ [0 for i in xrange(self.total_orb)] for j in xrange(self.total_orb)] # defining the dimension of the matrix.
        scc_energy_init = 0.0  # initial scc energy
        tmp_energy_init = 0.0  # initial tmp_energy
        tmp_energy_final = 0.0  # final tmp_energy      
        #tmp_scc_dftb_total_energy_init = 0.0
        #tmp_scc_dftb_total_energy_final = 0.0
        net_charge_init  = [0.0 for i in xrange(self.n_atom)]
        net_charge_Final = [0.0 for i in xrange(self.n_atom)]
        
#        diff = 0.0
        scc_count = 0
        scc_count_final = 0

        print ""
        print "--- Beginning of SCC-DFTB information ---"
        print ""

        while (scc_count >= 0):

              scc_diagonalization = Diagonalization(self.scc_dftb_matrix, 
                                                    self.S_dftb,
                                                    self.total_orb)
              scc_diagonalization.dftb_matrix_diagonalization()

              eigen_vectors_final = scc_diagonalization.eigen_vectors
              eigen_values_final =  scc_diagonalization.eigen_values

              density =  DensityMatrix(eigen_vectors_final,
                                       eigen_values_final,
                                       self.total_orb,
                                       self.n_electrons)

              density.build_P_matrix_from_C()
              self.P_dftb_final = density.P_dftb
              self.P_dftb_e_final = density.P_dftb_e

              analysis = Analysis(self.atom_type,
                                  self.xyz_atom_symbols,
                                  self.max_angular_momentum,
                                  self.total_orb,
                                  self.n_electrons_atom,
                                  self.S_dftb,
                                  density.P_dftb)
              analysis.mulliken_analysis()

              self.net_charge = analysis.net_charge

              # Mixing of net charges
              if self.isMixing == True and scc_count > 0:
                  
                 self.net_charge = self.get_net_charge_mixing(self.net_charge, net_charge_init, self.scc_mixing_value)

              scc_matrix = SCCMatrix(self.slako_read_list,
                                     self.atom_type,
                                     self.atoms_type_combinations_list,
                                     self.n_atom,
                                     self.xyz_atom_symbols,
                                     self.dist_matrix,
                                     self.total_orb,
                                     self.atom_numb,
                                     self.H_dftb,
                                     self.S_dftb,
                                     self.net_charge)
              scc_matrix.get_hubbard_parameter()
#              scc_matrix.build_gamma_matrix()
              scc_matrix.build_gamma_matrix_from_C()
#              scc_matrix.build_scc_matrix()
              scc_matrix.build_scc_matrix_from_C()
              scc_matrix.build_scc_dftb_matrix()

              self.scc_matrix_without_S_final = scc_matrix.scc_matrix_without_S

              # Evaluating the electronic energy

              n_ocupied_orbitals = int(self.n_electrons/2)
              n_unocupied_orbitals = len(eigen_values_final) - n_ocupied_orbitals
              ocupied_orbitals_energies = []
        
              # Evaluating the scc energy
              scc_energy_final = 0.0
              scc_dftb_total_energy_final = 0.0

              scc_energy_final = self.get_scc_energy(self.gamma, self.net_charge)

              self.get_H0_energy(self.H_dftb, density.P_dftb)
              scc_dftb_total_energy_final = self.H0_energy + scc_energy_final + self.repulsion_energy

              print "-----------------------------------"
              print "SCC step:", scc_count + 1
              print "-----------------------------------"

              if self.debug == True:
                  print ""
                  print "Occupied orbital energies:"
                  count = 0
                  for i in xrange(int(self.n_electrons/2)):
                      count = count + 1
                      if count < 5:
                          print ("% .8f" % (eigen_values_final[i])),
                      else:
                          count = 0
                          print ("% .8f" % (eigen_values_final[i]))
                  print ""
                  print ""
              
              print "Energies (hartree):"
              print ("H0               : % .12f" % (self.H0_energy))
              print ("Escc             : % .12f" % (scc_energy_final))
              print ("H0 + Escc        : % .12f" % (self.H0_energy + scc_energy_final))
              print ("H0 + Escc + Erep : % .12f" % (scc_dftb_total_energy_final))
              print ""

              self.scc_dftb_matrix = scc_matrix.scc_dftb_matrix
#              self.net_charge = analysis.net_charge

              #tmp_scc_dftb_total_energy_final = scc_dftb_total_energy_final
              net_charge_final = self.net_charge

#              diff = abs(tmp_scc_dftb_total_energy_final - tmp_scc_dftb_total_energy_init)
              net_charge_rms = self.get_net_charge_rms(net_charge_final, net_charge_init)

              print ("Net charges RMS  : % .8e" % (net_charge_rms))
              print ""

#              if diff <= self.scc_tolerance:
              if net_charge_rms <= self.scc_tolerance:
                 scc_count_final = scc_count + 1
                 scc_count = -1
                 self.scc_energy = scc_energy_final
                 
                 self.scc_dftb_total_energy = scc_dftb_total_energy_final
                 

              else:
                 #tmp_scc_dftb_total_energy_init = tmp_scc_dftb_total_energy_final         
                 net_charge_init = net_charge_final
                 scc_count += 1

              if scc_count > self.scc_step:

                 print ""
                 print "  The program stopped in the", scc_count, "scc cycle"
                 print "  SCC did not converge with the desired scc tolerance,", self.scc_tolerance, ", and steps,", self.scc_step
                 print ""

                 sys.exit()
 
        print "-------------------------------------------"
        print "  SCC converged in", scc_count_final, "iterations :) "
        print "-------------------------------------------"
        print ""

        self.get_H0_energy(self.H_dftb, density.P_dftb)
        
        self.scc_dftb_total_energy = self.H0_energy + self.scc_energy + self.repulsion_energy

        print ("  Band energy (non-SCC) (hartree) = % .12f" % (self.band_energy))
        print ("  H0 energy             (hartree) = % .12f" % (self.H0_energy))
        print ("  SCC energy            (hartree) = % .12f" % (self.scc_energy))
        print ("  Repulsion energy      (hartree) = % .12f" % (self.repulsion_energy))
        print ("  SCC-DFTB total energy (hartree) = % .12f" % (self.scc_dftb_total_energy))

        print ""
        print "--- End of SCC-DFTB information ---"
        print ""

    def get_scc_energy(self, gamma, net_charge):

        scc_energy = 0.0

        for in_atom in xrange(self.n_atom):
            for jn_atom in xrange(self.n_atom):

                scc_energy = scc_energy + 0.5*gamma[in_atom][jn_atom]*net_charge[in_atom]*net_charge[jn_atom]

        return scc_energy

    def get_H0_energy(self, H_dftb, P_dftb):
        
        tmp = np.dot(P_dftb, H_dftb)
        tmp2 = 0.0
        for i in xrange(self.total_orb):
            tmp2 = tmp2 + tmp[i][i]
        
        self.H0_energy = tmp2
        return self.H0_energy
        
    def get_net_charge_rms(self, net_charge, net_charge_old):
                
        #print "net_charge"
        #print net_charge
        #print net_charge_old
        
        net_charge_rms = 0.0
        average = 0.0
        for i in xrange(self.n_atom):
            average = average + (net_charge[i] - net_charge_old[i])**2.0

        net_charge_rms = math.sqrt((1.0/self.n_atom)*average)
        
        #print "net_charge_rms", net_charge_rms
        
        return net_charge_rms

    def get_net_charge_mixing(self, net_charge, net_charge_old, scc_mixing_value):
        
        for i in xrange(self.n_atom):
            
            net_charge[i] = net_charge[i]*(1.0 - scc_mixing_value) + net_charge_old[i]*(scc_mixing_value)

        return net_charge

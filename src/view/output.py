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

class OutPut(object):

    def __init__(self, passed_list_of_input_keywords,
                 type_calculation, 
                 path_slako,
                 atom_type,
                 max_angular_momentum,
                 max_all_momentum,
                 n_atom,
                 total_charge,
                 n_electrons,
                 scc_tolerance,
                 scc_step,
                 scc_mixing_value,
                 xyz_atom_symbols,
                 xyz,
                 isGeometryOptimization,
                 geometry_optimization_method,
                 geometry_optimization_tolerance,
                 geometry_optimization_steps,
                 geometry_out_file,
                 debug):
  
        self.passed_list_of_input_keywords = passed_list_of_input_keywords
        self.type_calculation = type_calculation
        self.path_slako = path_slako
        self.atom_type = atom_type
        self.max_angular_momentum = max_angular_momentum
        self.max_all_momentum = max_all_momentum
        self.n_atom =  n_atom
        self.total_charge = total_charge
        self.n_electrons  = n_electrons
        self.scc_tolerance = scc_tolerance
        self.scc_step = scc_step
        self.scc_mixing_value = scc_mixing_value
        self.xyz_atom_symbols = xyz_atom_symbols
        self.xyz = xyz # xyz MUST net in bohr
        self.isGeometryOptimization = isGeometryOptimization
        self.geometry_optimization_method = geometry_optimization_method
        self.geometry_optimization_tolerance = geometry_optimization_tolerance
        self.geometry_optimization_steps = geometry_optimization_steps
        self.geometry_out_file = geometry_out_file 
        self.debug = debug
        
        self.factor_bohr_to_angstrom = 1.0/1.889725989
 
    def print_output_header(self):

        print ""
        print "-----------------------------------------------------------------------------"
        print "-----------------------------------------------------------------------------"
        print " DFTB-Py is a DFTB/SCC-DFTB electronic structure code written in Python and C"
        print "-----------------------------------------------------------------------------"
        print "-----------------------------------------------------------------------------\n"

        print "Coder: Dr. Maicon Pierre Lourenco (MPL)."
        print "Collaborator: Mauricio Chagas da Silva."
 
        print "This code can be downloaded, distributed and modified for FREE under the GPL (GNU General Public License) license terms."
        print "Anyway, seize your day and learn a bit of the DFTB/SCC-DFTB method exploring the nice features of Python/C."
        print ""
        print "-------------------------------------------------------------------------"
        print "----------------------------- BEGIN OF OUTPUT ---------------------------"
        print "-------------------------------------------------------------------------"
        print ""
        print "--- Beginning of input information ---"
        print ""

        print "  List of valid keywords read in the input:"
        print " ", self.passed_list_of_input_keywords
        print ""
        print "  Type of calculation:", self.type_calculation
        if self.isGeometryOptimization == True:
            print "  Geometry optimization method:", self.geometry_optimization_method
            print "  Optimization tolerance:", self.geometry_optimization_tolerance, "(default is 1.0E-4)"
            print "  N. optimization steps:", self.geometry_optimization_steps, "(default is 500)"
        else:
            print "  Single point"
        print "  Slako's path:", "'" + self.path_slako + "'"
        print "  Type of atoms:", self.atom_type
        print "  Maximum angular momentum for each atom type:", self.max_angular_momentum
        print "  Maximum angular momentum of all atoms type:", self.max_all_momentum
        print "  Number of atoms:", int(self.n_atom)
        print "  Total charge of the system:", self.total_charge
        print "  Total number of valence electrons:", int(self.n_electrons)
        if self.type_calculation == "SCC-DFTB":
            print "  SCC tolerance:", self.scc_tolerance, "(default is 1.0E-8)"
            print "  SCC maximum number of steps:", self.scc_step, "(default is 100)"
            print "  Net charge mixing:", self.scc_mixing_value, "(default is 0.2)"

        print ""
        print "  Initial cartesian coordinates in bohr"
        print ""

        for row in xrange(self.n_atom):
            print('%6s % .8f % .8f % .8f' % (self.xyz_atom_symbols[row], \
                                             self.xyz[row,0], \
                                             self.xyz[row,1], \
                                             self.xyz[row,2]))
        print ""

        print "  Initial cartesian coordinates in angstrom"
        print ""

        for row in xrange(self.n_atom):
            print('%6s  % .8f % .8f % .8f' % (self.xyz_atom_symbols[row], \
                                             self.xyz[row,0]*self.factor_bohr_to_angstrom, \
                                             self.xyz[row,1]*self.factor_bohr_to_angstrom, \
                                             self.xyz[row,2]*self.factor_bohr_to_angstrom))
        print ""
        print "--- End of input information ---"

    def print_DFTB_data(self, N_Mulliken,
                         net_charge,
                         eigen_values,
                         band_energy,
                         repulsion_energy,
                         dftb_total_energy,
                         geometryOptimization_step):

        print ""
        print "--- Beginning of standard DFTB information ---"
        print ""

        if geometryOptimization_step == 0:

            if self.debug == True:
    
               print "  DFTB Mulliken population analysis of atoms "
               print ""
    
               for in_atom in xrange(self.n_atom):
                   print ("%6s % .8f" % (self.xyz_atom_symbols[in_atom], N_Mulliken[in_atom]))
               print ""
    
            print "  DFTB net charge from Mulliken population analysis of atoms "
            print ""
    
            for in_atom in xrange(self.n_atom):
                print ("%6s % .8f" % (self.xyz_atom_symbols[in_atom], net_charge[in_atom]))
    
            print ""

        if self.debug == True:
            print "  Occupied orbital energies:"
            count = 0
            for i in xrange(int(self.n_electrons/2)):
                count = count + 1
                if count < 5:
                     print ("  % .8f" % (eigen_values[i])),
                else:
                     count = 0
                     print ("  % .8f" % (eigen_values[i]))
            print ""
            print ""

        print ("  DFTB electronic energy (hartree) = % .12f" % (band_energy))
        print ("  Repulsion energy       (hartree) = % .12f" % (repulsion_energy))
        print ("  DFTB total energy      (hartree) = % .12f" % (dftb_total_energy))
        print ""

        print "--- End of standard DFTB information ---\n"

    def print_forces_data(self, total_forces, type_calculation):
       
        if type_calculation == "DFTB":
            print "Electronic + Repulsion forces:"
        else:
            print "Electronic + SCC + Repulsion forces:"
        
        for iatom in xrange(self.n_atom):
            print ("%2s % .12f % .12f % .12f, module:% .12f" % (self.xyz_atom_symbols[iatom], \
                                                                total_forces[iatom][0], \
                                                                total_forces[iatom][1], \
                                                                total_forces[iatom][2], \
                                                                sqrt(total_forces[iatom][0]**2 + total_forces[iatom][1]**2 + total_forces[iatom][2]**2)))

    def print_geometry_optimization_data_on_file(self, new_xyz,
                                    net_charge,
                                    geometryOptimization_step):
       
        self.geometry_out_file.write("%i" % self.n_atom + "\n")
        self.geometry_out_file.write("%s %i" % ("Coordinates (ang.) and the net_charge of the atoms. Step" , geometryOptimization_step ) + "\n")
        for i in xrange(self.n_atom):
            self.geometry_out_file.write("%2s % .12f % .12f % .12f % .12f" % (self.xyz_atom_symbols[i], \
                                                                         new_xyz[i][0]*self.factor_bohr_to_angstrom, \
                                                                         new_xyz[i][1]*self.factor_bohr_to_angstrom, \
                                                                         new_xyz[i][2]*self.factor_bohr_to_angstrom,
                                                                         net_charge[i]) + "\n")
    
    def print_geometry_optimized_data(self, new_xyz,
                                            net_charge,
                                            total_forces,
                                            geometryOptimization_step,
                                            isGeometryConverged):
       
        print "New coordinates got from the previous forces:"
        print ("%26s %36s" % ("angstrom", "bohr"))
       
        for i in xrange(self.n_atom):
            print ("%2s % .8f  % .8f  % .8f    % .8f  % .8f  % .8f" % (self.xyz_atom_symbols[i], \
                                                                        new_xyz[i][0]*self.factor_bohr_to_angstrom, \
                                                                        new_xyz[i][1]*self.factor_bohr_to_angstrom, \
                                                                        new_xyz[i][2]*self.factor_bohr_to_angstrom, \
                                                                        new_xyz[i][0], new_xyz[i][1], new_xyz[i][2]))
       
        if isGeometryConverged == True:
                                      
            print ""
            print "Converged geometry (in angstrom):"
            for iatom in xrange(self.n_atom):
                print ("%2s % .8f % .8f % .8f" % (self.xyz_atom_symbols[i], \
                      new_xyz[iatom][0]*self.factor_bohr_to_angstrom, \
                      new_xyz[iatom][1]*self.factor_bohr_to_angstrom, \
                      new_xyz[iatom][2]*self.factor_bohr_to_angstrom))
           
            print ""
            print "Final forces on atoms of the converged geometry:"
            for iatom in xrange(self.n_atom):
                print ("%2s % .12f % .12f % .12f, module:% .12f" % (self.xyz_atom_symbols[iatom], \
                                                                    total_forces[iatom][0], \
                                                                    total_forces[iatom][1], \
                                                                    total_forces[iatom][2], \
                                                                    sqrt(total_forces[iatom][0]**2 + total_forces[iatom][1]**2 + total_forces[iatom][2]**2)))

               
    def print_final_electronic_data(self, net_charge):
    
        print "  Net charge from Mulliken population analysis of atoms "
        print ""
    
        for in_atom in xrange(self.n_atom):
            print ("%6s % .8f" % (self.xyz_atom_symbols[in_atom], net_charge[in_atom]))
        
        print ""
    
    def print_output_tail(self):

        print "-------------------------------------------------------------------------"
        print "------------------------------ END OF OUTPUT ----------------------------"        
        print "-------------------------------------------------------------------------\n"

        print "  We who cut mere stones must always be envisioning cathedrals. Quarry worker's creed."
        print ""
              

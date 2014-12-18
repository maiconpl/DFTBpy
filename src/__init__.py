#!/usr/bin/python
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

__version__ = "1.0.55"

from view.parser import InputReader
from view.parser import SlakoReader
from util.dftb_util import DFTBMatrixUtil
from model.dftb_matrix import DFTBMatrix
from model.diag import Diagonalization
from model.density import DensityMatrix
from model.analysis import Analysis
from model.scc_matrix import SCCMatrix
from model.scc_procedure import SCCProcedure
from model.energy import Energy
from view.output import OutPut
from model.forces import ForcesOnAtom
from model.geometry_optimization import GeometryOptimization
import sys

# Measuring the execution time:
import datetime
from numpy import matrixlib
# Initial time:
t0 = datetime.datetime.now()

#--------------------------------------------------------
# Creating the instances and the objects
#--------------------------------------------------------

if __name__ == "__main__":
   # Creating the object "InputReader" and calling its the methods. The order of the calls does matter.
   input_read = InputReader("inp.dftb")
   input_read.get_type_of_calculations()
   input_read.get_slako_path()
   input_read.get_atom_type()
   input_read.get_max_angular_momentum()
   input_read.get_max_all_momentum()
   input_read.get_number_of_atoms()
   input_read.get_total_charge()
   input_read.get_xyz()
   
   isGeometryConverged = True

   def get_slako_read_list(atom_type, max_angular_momentum, path_slako):

       n_atom_type = len(atom_type)
       slako_read_list = [ [0 for i in xrange(n_atom_type)] for j in xrange(n_atom_type)]
       atoms_type_combinations_list = [ [0 for i in xrange(n_atom_type)] for j in xrange(n_atom_type)]

       # Creating the instances to read the slako file for all atom pair combinations.
       # That's is used into others methods of this class or, my be, in other classes.
       for iatom_type in xrange(n_atom_type):
           for jatom_type in xrange(n_atom_type):
               slako_read_list[iatom_type][jatom_type] = SlakoReader(path_slako,
                                                                     atom_type[iatom_type],
                                                                     atom_type[jatom_type])
               slako_read_list[iatom_type][jatom_type].get_slako()
               
               # Beginning of handling the max_all_momentum of each atom type
               if iatom_type == jatom_type:
                   
                   zero_counter = 0
                   
                   for i in xrange(3):
    
                       if slako_read_list[iatom_type][iatom_type].orb_energy[i] == 0.0:

                           zero_counter = zero_counter + 1
                   
                   try:
                       if zero_counter == 3:
                           raise
                   
                   except:
                        
                        print ""
                        print "ERROR!!!"
                        print "Error in the SLAKO pair " + "'" + atom_type[iatom_type] + "-" + atom_type[iatom_type] + "'."
                        print "The valence orbital energies read in the SLAKO '" + atom_type[iatom_type] + "-" + atom_type[iatom_type] + "'"  " are: " + str(slako_read_list[iatom_type][iatom_type].orb_energy) + "."
                        print "They can't be all zero."
                        print "Please, check the correctness of the SLAKO."
                        sys.exit()
                        
                   try:
                       
                       if zero_counter == 2 and (max_angular_momentum[iatom_type] == "d" or max_angular_momentum[iatom_type] == "p"):                            
                            raise

                       if zero_counter == 1 and (max_angular_momentum[iatom_type] == "d"):
                            raise
                   
                   except:
                        
                        print 
                        print "ERROR!!!!"
                        print "The maximum angular momentum declared of the atom " \
                               + "'" + atom_type[iatom_type] + "' is " + "'" + max_angular_momentum[iatom_type] + "'."
                        print "But, the valence orbital energies read in the SLAKO pair '" + atom_type[iatom_type] + "-" + atom_type[iatom_type] + "'"  " are: " + str(slako_read_list[iatom_type][iatom_type].orb_energy) + "."
                        print "For the maximum angular momentum " + "'" +  max_angular_momentum[iatom_type] + "'", "the number of zero(s) can't be" + " '" + str(zero_counter) + "'."

                        if max_angular_momentum[iatom_type] == "p":
                            print "It ought to have just ONE orbital energy (the d one) with value '0.0'."
                        
                        if max_angular_momentum[iatom_type] == "d":
                            print "ALL the valence orbital energies ought to be DIFFERENT of 0.0."

                        sys.exit()
               # End of handling the max_all_momentum of each atom type
               
               atoms_type_combinations_list[iatom_type][jatom_type] = atom_type[iatom_type] + atom_type[jatom_type]

       return slako_read_list, atoms_type_combinations_list
   
   slako_read_list, atoms_type_combinations_list = get_slako_read_list(input_read.atom_type,
                                                                       input_read.max_angular_momentum, 
                                                                       input_read.path_slako)
   
   # Creating the objects related with the class "DFTBMatrix"
   matrixUtil = DFTBMatrixUtil(slako_read_list,
                       input_read.atom_type,
                       input_read.n_atom,
                       input_read.total_charge,
                       input_read.max_angular_momentum,
                       input_read.max_all_momentum,
                       input_read.xyz_atom_symbols,
                       input_read.xyz,
                       input_read.debug)

   matrixUtil.build_data_lists()
   matrixUtil.build_orbitals_index_range_list_per_atom()

   H_list = None
   S_list = None
   dist_list = None
   H_list, S_list, dist_list, dist_step_list = matrixUtil.build_H_S_dist_list()
   
#    fit_H_list = None
#    fit_S_list = None
#    H_cutoff_list = None
#    fit_H_list, fit_S_list = matrixUtil.build_skf_fit_list(dist_list,
#                                                           H_list,
#                                                           S_list)

#    To print the fit of the matrix elements
#   matrixUtil.print_skf_fit(dist_list, H_list, S_list, fit_H_list, fit_S_list)

   H_cutoff_list = matrixUtil.build_integral_cutoff_list(H_list, dist_list)
   
#    matrixUtil.test_polint_fit_(dist_list[0][0][9],
#                                dist_step_list[0][0],
#                                H_list[0][0][9],
#                                S_list[0][0][9])
   
   matrixUtil.get_number_of_electrons()

   xyz = input_read.xyz
   xyz_ang = input_read.xyz_ang
   dist_matrix, dist_x, dist_y, dist_z = matrixUtil.build_dist_matrix(input_read.n_atom, xyz)

   geometry_out_file = open("geom.xyz", "w")
   output = OutPut(input_read.passed_list_of_input_keywords,
                   input_read.type_calculation,
                   input_read.path_slako,
                   input_read.atom_type,
                   input_read.max_angular_momentum, 
                   input_read.max_all_momentum,
                   input_read.n_atom,
                   input_read.total_charge,
                   matrixUtil.n_electrons,
                   input_read.scc_tolerance,
                   input_read.scc_step,
                   input_read.scc_mixing_value,
                   input_read.xyz_atom_symbols,
                   xyz,
                   input_read.isGeometryOptimization,
                   input_read.geometry_optimization_method,
                   input_read.geometry_optimization_tolerance,
                   input_read.geometry_optimization_steps,
                   geometry_out_file,
                   input_read.debug)
   output.print_output_header()

   count_opt_step = 0
   isGeometryConverged = False
   xyz_old = None
   total_forces_old = None
   total_energy = 0.0
   total_energy_old = 0.0

   while True: # concerning the geometry optimization steps

       matrix = DFTBMatrix(dist_x,
                           dist_y,
                           dist_z,
                           dist_matrix,
                           matrixUtil.atom_symb,
                           matrixUtil.atom_numb,
                           matrixUtil.atom_type_numb,
                           matrixUtil.max_all_momentum,
                           matrixUtil.atom_orb_sym,
                           matrixUtil.total_orb,
                           matrixUtil.H_diag_dftb,
                           dist_list,
                           dist_step_list,
                           H_list, 
                           S_list,
                           H_cutoff_list,
                           input_read.debug)
       
       matrix.build_H_S_matrix()
       if input_read.debug == True:
           matrix.print_H_S_matrix()
    
       diagonalization = Diagonalization(matrix.H_dftb,
                                         matrix.S_dftb, 
                                         matrix.total_orb)
       diagonalization.dftb_matrix_diagonalization()
    
       density = DensityMatrix(diagonalization.eigen_vectors, 
                               diagonalization.eigen_values,
                               matrix.total_orb,
                               matrixUtil.n_electrons)
       density.build_P_matrix_from_C()
    
       analysis = Analysis(input_read.atom_type,
                           input_read.xyz_atom_symbols,
                           input_read.max_angular_momentum,
                           matrix.total_orb,
                           matrixUtil.n_electrons_atom,
                           matrix.S_dftb,
                           density.P_dftb)
       analysis.mulliken_analysis()
    
       energy = Energy(slako_read_list,
                       matrixUtil.atom_type,
                       atoms_type_combinations_list,
                       diagonalization.eigen_values,
                       matrixUtil.n_electrons,
                       input_read.n_atom,
                       input_read.xyz_atom_symbols,
                       dist_matrix)
       energy.get_band_energy()
       energy.get_repulsion_energy()
       energy.get_dftb_total_energy()
       
       if input_read.type_calculation == "DFTB" \
          or (input_read.type_calculation == "SCC-DFTB" and count_opt_step == 0):
           output.print_DFTB_data(analysis.N_Mulliken,
                        analysis.net_charge,
                        diagonalization.eigen_values,
                        energy.band_energy,
                        energy.repulsion_energy,
                        energy.get_dftb_total_energy(),
                        count_opt_step)

       if input_read.type_calculation == "DFTB":
           net_charge = analysis.net_charge
           output.print_geometry_optimization_data_on_file(xyz, net_charge, count_opt_step)

       if input_read.type_calculation == "DFTB" and input_read.isGeometryOptimization == False:
          break
    
       if input_read.type_calculation == "DFTB" and isGeometryConverged == True:
           output.print_final_electronic_data(analysis.net_charge)
           break
       
       if input_read.type_calculation == "DFTB" and count_opt_step >= input_read.geometry_optimization_steps:

           print "WARNING! Geometry did not converge in " + str(count_opt_step) + " iterations.\n"
           break
       
       if input_read.type_calculation == "DFTB" and input_read.isGeometryOptimization == True:
           
           total_energy = energy.dftb_total_energy
           print "DFTB forces on atoms:"
           
           forces = ForcesOnAtom(dist_list,
                                 dist_step_list,
                                 H_list, 
                                 S_list,
                                 slako_read_list,
                                 dist_x,
                                 dist_y,
                                 dist_z,
                                 dist_matrix,
                                 input_read.n_atom,
                                 input_read.xyz_atom_symbols,
                                 input_read.atom_type,
                                 matrixUtil.atom_symb,
                                 matrixUtil.atom_numb,
                                 matrixUtil.atom_type_numb,
                                 atoms_type_combinations_list,
                                 matrixUtil.orbitals_index_range_list_per_atom,
                                 matrixUtil.max_all_momentum,
                                 matrixUtil.atom_orb_sym,
                                 matrixUtil.total_orb,
                                 density.P_dftb,
                                 density.P_dftb_e,
                                 analysis.net_charge,
                                 None,
                                 False,
                                 input_read.debug)     
           forces.get_dftb_forces()
    
       if input_read.type_calculation == "SCC-DFTB":
    
          scc_matrix = SCCMatrix(slako_read_list,
                                 matrixUtil.atom_type,
                                 atoms_type_combinations_list,                             
                                 input_read.n_atom,
                                 input_read.xyz_atom_symbols,
                                 dist_matrix,
                                 matrix.total_orb,
                                 matrix.atom_numb,
                                 matrix.H_dftb,
                                 matrix.S_dftb,
                                 analysis.net_charge)
          scc_matrix.get_hubbard_parameter()
          scc_matrix.build_gamma_matrix_from_C()
          scc_matrix.build_scc_matrix_from_C()
          scc_matrix.build_scc_dftb_matrix()
    
          scc_procedure = SCCProcedure(slako_read_list,
                                       atoms_type_combinations_list,
                                       input_read.scc_step,
                                       input_read.scc_tolerance,
                                       input_read.atom_type,
                                       input_read.xyz_atom_symbols,
                                       dist_matrix,
                                       input_read.max_angular_momentum, 
                                       matrixUtil.n_electrons_atom,
                                       scc_matrix.scc_dftb_matrix,
                                       matrix.H_dftb,
                                       matrix.S_dftb,
                                       scc_matrix.gamma,
                                       matrix.total_orb,
                                       matrix.atom_numb,
                                       input_read.n_atom, 
                                       matrixUtil.n_electrons,
                                       energy.get_band_energy(),
                                       energy.repulsion_energy,
                                       input_read.scc_mixing_value,
                                       input_read.debug)
          scc_procedure.scc_procedure()

       if input_read.type_calculation == "SCC-DFTB":
           net_charge = scc_procedure.net_charge
           output.print_geometry_optimization_data_on_file(xyz, net_charge, count_opt_step)
    
       if input_read.type_calculation == "SCC-DFTB" and input_read.isGeometryOptimization == False:
           break

       if input_read.type_calculation == "SCC-DFTB" and isGeometryConverged == True:
           output.print_final_electronic_data(scc_procedure.net_charge)
           break

       if input_read.type_calculation == "SCC-DFTB" and count_opt_step >= input_read.geometry_optimization_steps:

           print "WARNING! Geometry did not converge in " + str(count_opt_step) + " iterations.\n"
           break

       if input_read.type_calculation == "SCC-DFTB" and input_read.isGeometryOptimization == True:

          total_energy = scc_procedure.scc_dftb_total_energy     
          print "SCC-DFTB forces on atoms:"

          forces = ForcesOnAtom(dist_list,
                                dist_step_list,
                                H_list,
                                S_list,
                                slako_read_list,
                                dist_x,
                                dist_y,
                                dist_z,
                                dist_matrix,
                                input_read.n_atom,
                                input_read.xyz_atom_symbols,
                                input_read.atom_type,
                                matrixUtil.atom_symb,
                                matrixUtil.atom_numb,
                                matrixUtil.atom_type_numb,
                                atoms_type_combinations_list,
                                matrixUtil.orbitals_index_range_list_per_atom,
                                matrixUtil.max_all_momentum,
                                matrixUtil.atom_orb_sym,
                                matrixUtil.total_orb,
                                scc_procedure.P_dftb_final,
                                scc_procedure.P_dftb_e_final,
                                scc_procedure.net_charge,
                                scc_procedure.scc_matrix_without_S_final,
                                True,
                                input_read.debug)
          forces.get_dftb_forces()
          
       output.print_forces_data(forces.total_forces,
                                input_read.type_calculation)

       print ""
       print "-------------------------------------------"
       print "Geometry optimization step:", count_opt_step + 1
       print "-------------------------------------------\n"

       optimization = GeometryOptimization(input_read.n_atom,
                                           forces.total_forces,
                                           input_read.geometry_optimization_tolerance)

       if input_read.geometry_optimization_method == "ssd":
    
           weight = [ [1.0 for i in xrange(3)] for j in xrange(input_read.n_atom)]
           optimization.do_simple_steepest_descent(xyz, weight, total_energy, total_energy_old)
       
       if input_read.geometry_optimization_method == "sd":

           if count_opt_step == 0:
               weight = [ [0.5 for i in xrange(3)] for j in xrange(input_read.n_atom)]
               weight_old = [ [0.5 for i in xrange(3)] for j in xrange(input_read.n_atom)]

           weight = optimization.do_steepest_descent(xyz,
                                            xyz_old,
                                            forces.total_forces,
                                            total_forces_old,
                                            weight,
                                            weight_old,
                                            total_energy,
                                            total_energy_old,
                                            count_opt_step)
           weight_old = weight
       isGeometryConverged = optimization.isGeometryConverged

       # updating the new coordinates, forces and energies. That's a crucial part of the program.
       xyz_old = xyz
       total_forces_old = forces.total_forces
       xyz = optimization.xyz_new
       total_energy_old = total_energy
       dist_matrix, dist_x, dist_y, dist_z = matrixUtil.build_dist_matrix(input_read.n_atom, xyz)

       count_opt_step =  count_opt_step + 1
       
       if count_opt_step > 0:
           output.print_geometry_optimized_data(xyz,
                                                net_charge,
                                                forces.total_forces, 
                                                count_opt_step,
                                                isGeometryConverged)
         
   output.print_output_tail()
   geometry_out_file.close()
   
#-------------------------------------------
# Measuring the execution time
#-------------------------------------------

delta_t = datetime.datetime.now() - t0

#Final time:

print "  Execution time:", delta_t

# Date and time of execution
print ""
print "  Date and time of the execution:", t0.strftime("%Y-%m-%d, %H:%M")

# version
print ("")
print ("VERSION OF DFTBPy PROGRAM: DFTBPy-" + __version__ + " Released on 17th of December 2014.")

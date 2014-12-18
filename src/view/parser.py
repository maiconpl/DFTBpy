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
import sys
import os

class InputReader(object):

    def __init__(self, input_name):
 
        self.input_name = input_name
        self.lines = None
        self.type_calculation = None
        self.scc_step = None
        self.scc_tolerance = None
        self.path_slako = None
        self.atom_type = []
        self.max_angular_momentum = None
        self.max_all_momentum = None # max angular momentum of all atoms type
        self.n_atom = 0
        self.total_charge = 0
        self.cartesian_type = None
        self.xyz_atom_symbols = []
        self.xyz  = [] # in bohr
        self.xyz_ang  = [] # in angstrom
        self.isGeometryOptimization = False
        self.geometry_optimization_method = None
        self.geometry_optimization_steps = 0.0
        self.geometry_optimization_tolerance = 0.0
        self.scc_mixing_value = 0.2 # default value
        self.debug = None
        self.passed_list_of_input_keywords = []

    def read_input(self):

        try:
            fileIn = file(self.input_name, "r")
        except:
            
            print ""
            print "ERROR!!!"
            print "The input file" + " '" + self.input_name + "' " + "doesn't exist.\n" 
            
            sys.exit()

        self.lines = fileIn.readlines()

        fileIn.close()
         
        return self.lines

    def get_type_of_calculations(self):

        tmp = str(self.read_input()[2]).split()

        try:

            self.type_calculation = tmp[0]
            
            if len(self.type_calculation) == 0:
                raise

        except:
            print ""
            print "ERROR!!!"
            print "The type of calculation " + "'" + str(self.type_calculation) + "'" + " is not correct."
            print "The possible error is the blank line."
            print "Please, check if whole the input is set correctly."
            print "See the manual for instructions or some input template."

            sys.exit()
        
        self.passed_list_of_input_keywords.append(self.type_calculation)

        try:
        
            if self.type_calculation != "DFTB" and self.type_calculation != "SCC-DFTB":
                raise
        
        except:
                print ""
                print " ERROR!!!"
                print " This kind of calculation" + " '" + self.type_calculation + "' " + "isn't valid."
                print " It ought to be: DFTB OR SCC-DFTB."
                print ""

                sys.exit()

        if self.type_calculation == "SCC-DFTB":
            
           self.scc_step = 100 # default
           self.scc_tolerance = 10.0**(-8) # default
           self.scc_mixing_value = 0.2 # default

           try:
           
              for iString in xrange(len(tmp)):
    
                    if tmp[iString] == "scc_tol":
        
                        self.passed_list_of_input_keywords.append(tmp[iString])
        
                        if tmp[iString + 1] == "=":
        
                            self.scc_tolerance = float(tmp[iString + 2])
                        
                        elif tmp[iString + 1] != "=":
        
                            if isinstance(float(tmp[iString + 1].split("=")[1]), float) == True:
        
                                self.scc_tolerance = float(tmp[iString + 1].split("=")[1])
                            
                            else:
                        
                                self.scc_tolerance = float(tmp[iString + 1])
        
                    elif tmp[iString].strip() == "scc_tol=":
        
                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])
                        self.scc_tolerance = float(tmp[iString + 1])
        
                    elif tmp[iString].split("=")[0] == "scc_tol":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])         
                        self.scc_tolerance = float(tmp[iString].split("=")[1])
        
                    if tmp[iString] == "scc_step":

                        self.passed_list_of_input_keywords.append(tmp[iString])                 
                        
                        if tmp[iString + 1] == "=":
        
                            self.scc_step = float(tmp[iString + 2])
                            self.scc_step = int(self.scc_step)
                        
                        elif tmp[iString + 1] != "=":
        
                            if isinstance(float(tmp[iString + 1].split("=")[1]), float) == True:
        
                                self.scc_step = float(tmp[iString + 1].split("=")[1])
                                self.scc_step = int(self.scc_step)
                            
                            else:
                        
                                self.scc_step = float(tmp[iString + 1])
                                self.scc_step = int(self.scc_step)
        
                    elif tmp[iString].strip() == "scc_step=":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])        
                        self.scc_step = float(tmp[iString + 1])
                        self.scc_step = int(self.scc_step)
        
                    elif tmp[iString].split("=")[0] == "scc_step":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])         
                        self.scc_step = float(tmp[iString].split("=")[1])
                        self.scc_step = int(self.scc_step)

                    if tmp[iString] == "mix":

                        self.passed_list_of_input_keywords.append(tmp[iString])         
        
                        if tmp[iString + 1] == "=":
        
                            self.scc_mixing_value = float(tmp[iString + 2])
                        
                        elif tmp[iString + 1] != "=":
        
                            if isinstance(float(tmp[iString + 1].split("=")[1]), float) == True:
        
                                self.scc_mixing_value = float(tmp[iString + 1].split("=")[1])
                            
                            else:
                        
                                self.scc_mixing_value = float(tmp[iString + 1])
        
                    elif tmp[iString].strip() == "mix=":
        
                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])         
                        self.scc_mixing_value = float(tmp[iString + 1])
        
                    elif tmp[iString].split("=")[0] == "mix":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])         
                        self.scc_mixing_value = float(tmp[iString].split("=")[1])
                            
           except:
               
               print ""
               print "ERROR!!!"
               print "Error in" + " '" + self.input_name + "'."
               print "Please, check the general input, based on the manual."
               print "Also, you can check the third line of the input and see if the keywords are set like this:\n"
               print "SCC-DFTB scc_step=100.0 scc_tol=1.0E-8 mix=0.2 opt=sd geom_tol=1.0E-4 opt_step=100.0"
               print ""
               print "Probably, the error is/are in the keyword(s) involving the SCC."
               print "They ought to be set like this:"
               print ""
               print "scc_step=100.0 scc_tol=1.0E-8 mix=0.2"
               print ""
               
               sys.exit()

        index = 0

        if "opt" in tmp or "opt=" in tmp or "opt=sd" in tmp or "opt=ssd" in tmp or "opt=gc" in tmp:

            self.geometry_optimization_method = "sd" # default
            self.geometry_optimization_steps = 500 # default
            self.geometry_optimization_tolerance = 10.0**(-4) # default
            self.isGeometryOptimization = True

            self.passed_list_of_input_keywords.append("opt")

            try:
               
                for iString in xrange(len(tmp)):
    
                    if tmp[iString] == "opt":

                        self.geometry_optimization_method = "sd" # default                    
    
                    elif tmp[iString] == "opt=sd" or tmp[iString] == "opt=ssd" or tmp[iString] == "opt=gc":
    
                        self.geometry_optimization_method = tmp[iString].split("=")[1]
                    
                    elif tmp[iString] == "opt" or tmp[iString] == "opt=":
                        
                        if tmp[iString + 1] == "=":
                            
                            index = 2
                        
                        else:
            
                            index = 1
                            
                        if tmp[ tmp.index(tmp[iString]) + index].strip() == "sd" or \
                           tmp[ tmp.index(tmp[iString]) + index].strip() == "=sd":
                                
                            self.geometry_optimization_method = "sd"
                                
                        elif tmp[ tmp.index(tmp[iString]) + index].strip() == "ssd" or \
                             tmp[ tmp.index(tmp[iString]) + index].strip() == "=ssd":
    
                            self.geometry_optimization_method = "ssd"
    
                        elif tmp[ tmp.index(tmp[iString]) + index].strip() == "gc" or \
                             tmp[ tmp.index(tmp[iString]) + index].strip() == "=gc":
    
                            self.geometry_optimization_method = "gc"
                                
                        else:
                            print "ERROR!!!"
                            print "The 'opt' option " + "'" + tmp[ tmp.index(tmp[iString]) + index ] + "'" + " is not valid." 
                            print "The valid options is: steepest descent 'sd', simpler steepest descent 'ssd' and gradient conjugated 'gc'." 
                            sys.exit()
    
                    if tmp[iString] == "geom_tol":
    
                        self.passed_list_of_input_keywords.append(tmp[iString])
                        
                        if tmp[iString + 1] == "=":
    
                            self.geometry_optimization_tolerance = float(tmp[iString + 2])
                        
                        elif tmp[iString + 1] != "=":
    
                            if isinstance(float(tmp[iString + 1].split("=")[1]), float) == True:
    
                                self.geometry_optimization_tolerance = float(tmp[iString + 1].split("=")[1])
                            
                            else:
                        
                                self.geometry_optimization_tolerance = float(tmp[iString + 1])
    
                    elif tmp[iString].strip() == "geom_tol=":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])
                        self.geometry_optimization_tolerance = float(tmp[iString + 1])
    
                    elif tmp[iString].split("=")[0] == "geom_tol":
                        
                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])     
                        self.geometry_optimization_tolerance = float(tmp[iString].split("=")[1])
    
                    if tmp[iString] == "opt_step":

                        self.passed_list_of_input_keywords.append(tmp[iString])
    
                        if tmp[iString + 1] == "=":
    
                            self.geometry_optimization_steps = float(tmp[iString + 2])
                            self.geometry_optimization_steps = int(self.geometry_optimization_steps)
                        
                        elif tmp[iString + 1] != "=":
    
                            if isinstance(float(tmp[iString + 1].split("=")[1]), float) == True:
    
                                self.geometry_optimization_steps = float(tmp[iString + 1].split("=")[1])
                                self.geometry_optimization_steps = int(self.geometry_optimization_steps)
                            
                            else:
                        
                                self.geometry_optimization_steps = float(tmp[iString + 1])
                                self.geometry_optimization_steps = int(self.geometry_optimization_steps)
    
                    elif tmp[iString].strip() == "opt_step=":

                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])    
                        
                        self.geometry_optimization_steps = float(tmp[iString + 1])
                        self.geometry_optimization_steps = int(self.geometry_optimization_steps)
    
                    elif tmp[iString].split("=")[0] == "opt_step":
                        
                        self.passed_list_of_input_keywords.append(tmp[iString].split("=")[0])
                        
                        self.geometry_optimization_steps = float(tmp[iString].split("=")[1])
                        self.geometry_optimization_steps = int(self.geometry_optimization_steps)

            except:
               
                print ""
                print "ERROR!!!"
                print "Error in" + " '" + self.input_name + "'."
                print "Please, check the general input, based on the manual."
                print "Also, you can check the third line of the input and see if the keywords are set like this:\n"
                print "SCC-DFTB scc_step=100.0 scc_tol=1.0E-8 mix=0.2 opt=sd geom_tol=1.0E-4 opt_step=100.0"
                print ""
                print "Probably, the error is/are in the keyword(s) involving the geometry optimization."
                print "They ought to be set like this:"
                print ""
                print "opt=sd geom_tol=1.0E-8 opt_step=100.0"
                print ""
               
                sys.exit()

        if self.scc_mixing_value < 0.0 or self.scc_mixing_value > 1.0:
            
            print ""
            print "ERROR!!!"
            print "The mixing value is:", str(self.scc_mixing_value) + "."
            print "It ought to be between 0.0 and 1.0.\n"
            sys.exit()

    def get_slako_path(self):

        self.path_slako = str(self.read_input()[4]).strip()

        try:
            
            if os.path.isdir(self.path_slako) == False:
                raise 
            
        except:
            
            print ""
            print "ERROR!!!"
            print "The SLAKO's path" + " '" + self.path_slako + "' doesn't exist."
            print "Please, check if it was set correctly in the input.\n"
            
            sys.exit()

        return self.path_slako

    def get_atom_type(self):

        try:
            self.atom_type = self.read_input()[6].split()
            
            if len(self.atom_type) == 0:
                raise 

        except:
            print ""
            print "ERROR!!!"
            print "It was not possible to read the 'atoms type':" + str(self.atom_type) + "."
            print "Please, check if the input is set correctly."
            print "For instance, check the lines of the input."
            print "See the manual for instructions or some input template."

            sys.exit()

        return self.atom_type

    def get_max_angular_momentum(self):

        try:
            
            self.max_angular_momentum = self.read_input()[8].split()
            
            if len(self.max_angular_momentum) == 0:
                raise

        
        except:
            print ""
            print "ERROR!!!"
            print "It was not possible to read the 'maximum angular momentum' " + str(self.max_angular_momentum) + " of each atoms type."
            print "Please, check if the input is set correctly."
            print "For instance, check the lines of the input."
            print "See the manual for instructions or some input template."

            sys.exit()

        try:
            
            if len(self.max_angular_momentum) != len(self.atom_type):
                raise
        
        except:
            print ""
            print "ERROR!!!"
            print "The number of atoms type " + str(self.atom_type) + " IS NOT equal to the number of"
            print "the respective max angular momentum " + str(self.max_angular_momentum) + " of each atom."
            print "Please, check if the input is set correctly."
            print "For instance, also, check the lines of the input."
            print "See the manual for instructions or some input template."

            sys.exit()            

        return self.max_angular_momentum

    def get_max_all_momentum(self):

        count_s =0
        count_p =0
        count_d =0
        for imax_angular_momentum in self.max_angular_momentum:
            if imax_angular_momentum == "d":
               count_d += 1
            elif imax_angular_momentum == "p":
               count_p += 1
            elif imax_angular_momentum == "s":
               count_s += 1
            else:
                print ""
                print "ERROR!!!"
                print "The maximum angular momentum " + "'" + imax_angular_momentum + "' " + "is not supported or the symbol doesn't exist."
                print "It has to be just s, p or d valence orbitals.\n" 
                sys.exit()
        if count_d > 0:
           self.max_all_momentum = 2
        if count_d == 0 and count_p > 0:
           self.max_all_momentum = 1
        if count_d == 0 and count_p == 0 and count_s > 0:
           self.max_all_momentum = 0

        return self.max_all_momentum

    def get_number_of_atoms(self):

        try:
            self.n_atom = int(self.read_input()[10].split()[0])
        except:
            print ""
            print "ERROR!!!"
            print "It was not possible to read the number of atoms in the input."
            print "Please, check if the input is set correctly."
            print "For instance, check lines of the input."
            print "See the manual for instructions or some input template."
            sys.exit()

        # getting the debug variable

        if self.n_atom < 0:

           self.n_atom = abs(self.n_atom)

           self.debug = True

        return self.n_atom
    
    def get_total_charge(self):

        if len(self.read_input()[10].split()) == 2: # That works too
            
            self.total_charge = int(self.read_input()[10].split()[1])

    def get_cartesian_type(self): # bohr or angstrom
 
        words = self.read_input()[11].split()

        for word in words:
            if "angs" in word or "Angs" in word or "ANGS" in word:
               self.cartesian_type = "angstrom"
            if "bohr" in word or "Bohr" in word or "BOHR" in word:
               self.cartesian_type = "bohr"

        if self.cartesian_type != "angstrom" and self.cartesian_type != "bohr":
           
           print ""
           print "Error!!!"
           print "The 'angstrom' or 'bohr' units were not set correctly."
           print "Please, check the lines of the input."
           print "The xyz units -- 'angstrom' or 'bohr' -- have to be set like this:"
           print "# coordinates in angstrom"
           print "Si  -0.012513   -0.011043   -0.470230"
           print "Si   0.004534    0.004166    1.920228"
           sys.exit()

        return self.cartesian_type

    def get_xyz(self): # the xyz is read in angstrom

        self.get_cartesian_type()

        count_atoms = 0
        
        try:
            
            for line in self.lines[12:]:
                count_atoms = count_atoms + 1
                pieces = line.split()
                self.xyz_atom_symbols.append(pieces[0])
                self.xyz.append( tuple  ( (map( float, pieces[1:4]))))
            
        except:

            print ""            
            print "ERROR!!!"
            print "It was not possible to read the xyz coordinates in the input."
            print "The possible problems are:"
            print "(1) The number of lines and/or columns are not set correctly."
            print "(2) There is a blank line in the beginning, in the middle or, in most cases, in the end of the XYZ definition."
            print "(3) The is a string, instead of number, in the middle of the XYZ definition."
            sys.exit()

        self.xyz = array(self.xyz)

        try:
            if count_atoms != self.n_atom:
                raise
            
        except:
            print ""
            print "ERROR!!!"
            print "The number of defined atoms in the input " + "'" + str(self.n_atom) + "'"
            print "is not equal to the number of XYZ atomic coordinates: " + "'" + str(count_atoms) + "'."
            print "Please, fix it."
            sys.exit()

        
        tmp_list = []
        count_the_number_of_atoms_type_in_xyz = 0
        
        try:

            for iatom_type in xrange(len(self.atom_type)):            
                for j_atom in xrange(self.n_atom):
                    if self.atom_type[iatom_type] ==  self.xyz_atom_symbols[j_atom]:
                        count_the_number_of_atoms_type_in_xyz = count_the_number_of_atoms_type_in_xyz + 1
                        tmp_list.append(self.xyz_atom_symbols[j_atom])                        
                        
            if count_the_number_of_atoms_type_in_xyz != self.n_atom:
                
                raise
            
        except:

            print ""
            print "ERROR!!!"
            print "The defined atom types is: " + str(self.atom_type) + "."
            print "It is missing one or some atom type(s) present in the XYZ coordinate."
            print "You have to DEFINE it/them in the 'Atom Type' part of the input"
            print "or REMOVE the atomic coordinate(s) referent to the not defined atom type(s)."
            sys.exit()

        if self.cartesian_type == "angstrom":

           self.xyz_ang = self.xyz
           self.xyz_ang = array(self.xyz_ang)

           self.xyz = self.xyz_ang*1.889725989 # converting angstrom to bohr

        if self.cartesian_type == "bohr":

           self.xyz_ang = self.xyz/1.889725989 # converting bohr to angstrom

        return self.xyz # it must be returned always in bohr

class SlakoReader(object):
    """ This class reads the slako files"""

    def __init__(self, path, symbol_a, symbol_b):

        self.path = path
        self.symbol_a = symbol_a
        self.symbol_b = symbol_b

        self.dist_step = []  # Steps of the interatomic distances.
        self.n_step = []  # Number of steps of the interatomic distances into the slako file.
        self.orb_energy = []  # Orbital energy of the valence orbitals.
        self.hubbard_parameter = []  # Hubbard parameter of the d, p and s orbitals, respectively.
        self.orb_ocupation = []  # Ocupation of the d, p and s orbitals, respectively.
        self.atom_mass = []  # atomic mass.
        self.poly_coef = []  # not implemented yet
        self.cut_radius = []  # not implemented yet
        self.skf_H = []  # skf_H matrix which stores the integrals H integrals after the fitting of the slako elements
        self.skf_S = []  # skf_S matrix which stores the integrals S integrals after the fitting of the slako elements
 
    def get_slako(self):

        tmp = []
        skf_pair = self.symbol_a + '-' + self.symbol_b + ".skf"
        slako_address = self.path + skf_pair
        
        try:
            
            fileIn  = file(slako_address,"r")
        
        except:
            
            print ""
            print "ERROR!!!"
            print "It was not possible to open the file: " + "'" + slako_address + "'."
            print "Please, check if it exists or if it was set correctly in the input."
            print "The SLAKO's path ought to have '/' around, e.g:"
            print " /home/someone/SKF/"
            sys.exit()

        lines = fileIn.readlines()

        # reading the Homonuclear pairs
        if self.symbol_a == self.symbol_b:
           #
           tmp = lines[0].split()
           self.dist_step = float(tmp[0])
           self.n_step = int(tmp[1])
           #
           tmp = lines[1].split()
           self.orb_energy = (map (float, tmp[0:3]))
           self.hubbard_parameter = (map (float, tmp[4:7]))
           self.orb_ocupation = (map (float, tmp[8:10]))
           #
           tmp = lines[2].split()
           self.atom_mass = float(tmp[0])
           #
           self.poly_coef = (map (float, tmp[1:9]))
           self.cut_radius = float(tmp[9])
           #
           for line in lines[3:]:
               pieces = line.split() # split = separation
               self.skf_H.append(pieces[0:10])
               self.skf_S.append(pieces[10:20])

        # reading the Heteronuclear pairs
        else:
           tmp = lines[0].split()
           self.dist_step = float(tmp[0])
           self.n_step = int(tmp[1])
           #
           tmp = lines[1].split()
           self.poly_coef = (map (float, tmp[1:9]))
           self.cut_radius = float(tmp[9])
           #
           for line in lines[2:]:
               pieces = line.split() # split = separation
               self.skf_H.append(pieces[0:10])
               self.skf_S.append(pieces[10:20])
           #
 	fileIn.close()

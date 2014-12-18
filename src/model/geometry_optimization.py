'''
Created on Nov 12, 2014

@author: maicon
'''

class GeometryOptimization(object):

    def __init__(self, n_atom, total_forces, force_tolerance):
        
        self.n_atom = n_atom
        self.total_forces = total_forces
        self.force_tolerance = force_tolerance

        self.xyz_new = [ [0.0 for i in xrange(3)] for j in xrange(self.n_atom)]
        self.isGeometryConverged = False

    def do_simple_steepest_descent(self, xyz, weight, total_energy, total_energy_old):

        energy_diff = abs(total_energy - total_energy_old)
        
        count_converged_coordinates = 0
        
        for i in xrange(self.n_atom):
            
            if abs(self.total_forces[i][0]) >= self.force_tolerance:
                self.xyz_new[i][0] = xyz[i][0] - weight[i][0]*self.total_forces[i][0]
            else:
                count_converged_coordinates = count_converged_coordinates + 1
                self.xyz_new[i][0] = xyz[i][0]

            if abs(self.total_forces[i][1]) >= self.force_tolerance: 
                self.xyz_new[i][1] = xyz[i][1] -  weight[i][1]*self.total_forces[i][1]
            else:
                count_converged_coordinates = count_converged_coordinates + 1
                self.xyz_new[i][1] = xyz[i][1]

            if abs(self.total_forces[i][2]) >= self.force_tolerance:
                self.xyz_new[i][2] = xyz[i][2] -  weight[i][2]*self.total_forces[i][2]
            else:
                count_converged_coordinates = count_converged_coordinates + 1
                self.xyz_new[i][2] = xyz[i][2]

#        print "count_converged_coordinates", count_converged_coordinates
        if count_converged_coordinates == self.n_atom*3 and energy_diff <= 1.0E-6:
            self.isGeometryConverged = True
        
    def do_steepest_descent(self, xyz, xyz_old, 
                            total_forces, total_forces_old,
                            weight, weight_old,
                            total_energy, 
                            total_energy_old,
                            count_opt_step):
        
        energy_diff = abs(total_energy - total_energy_old)
        count_converged_coordinates = 0
        
        for i in xrange(self.n_atom):

            if count_opt_step > 0:

                if abs(total_forces[i][0]) >= self.force_tolerance:
                    weight[i][0] = abs(xyz[i][0] - xyz_old[i][0])/abs(total_forces[i][0] - total_forces_old[i][0])
                else:
                    weight[i][0] = weight_old[i][0]
                    count_converged_coordinates = count_converged_coordinates + 1

                if abs(total_forces[i][1]) >= self.force_tolerance:
                    weight[i][1] = abs(xyz[i][1] - xyz_old[i][1])/abs(total_forces[i][1] - total_forces_old[i][1])
                else:
                    weight[i][1] = weight_old[i][1]
                    count_converged_coordinates = count_converged_coordinates + 1

                if abs(total_forces[i][2]) >= self.force_tolerance:
                    weight[i][2] = abs(xyz[i][2] - xyz_old[i][2])/abs(total_forces[i][2] - total_forces_old[i][2])
                else:
                    weight[i][2] = weight_old[i][2]
                    count_converged_coordinates = count_converged_coordinates + 1

            self.xyz_new[i][0] = xyz[i][0] - weight[i][0]*total_forces[i][0]
            self.xyz_new[i][1] = xyz[i][1] - weight[i][1]*total_forces[i][1]
            self.xyz_new[i][2] = xyz[i][2] - weight[i][2]*total_forces[i][2]
            
#        print "count_converged_coordinates", count_converged_coordinates
        if count_converged_coordinates == self.n_atom*3 and energy_diff <= 1.0E-6:
            self.isGeometryConverged = True

        return weight
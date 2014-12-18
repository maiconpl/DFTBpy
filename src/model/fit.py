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

#from scipy.interpolate import interp1d, UnivariateSpline, InterpolatedUnivariateSpline, Rbf
from scipy.interpolate import InterpolatedUnivariateSpline
#from scipy.interpolate import interp1d
import sys

import os
_path = os.path.dirname(os.path.realpath(__file__)) + "/" + "../lib/polint.so"

from ctypes import cdll, byref, c_int, c_double, c_char
lib_polint = cdll.LoadLibrary(_path)

class Fit(object):

    def __init__(self):

        pass
   
    def slako_integrals_fit2(self, distance_grid, integral_grid):

        if len(integral_grid) == 0: # there is no value different form 0.0 for certain integral symmetry grid.

           matrix_element_fit = None

        if len(integral_grid) != 0:

           # Some tests indicate that InterpolatedUnivariateSpline and interp1d (kind="cubic") 
           # get the same results, which are better than interp1d(kind="linear").
           # But, InterpolatedUnivariateSpline method is much more faster (about 4 times)
           # than the interp1d (kind="cubic").
           
           matrix_element_fit = InterpolatedUnivariateSpline(distance_grid, integral_grid)
#            matrix_element_fit = interp1d(distance_grid, integral_grid, kind="cubic", bounds_error=False, fill_value=0.0)
#            matrix_element_fit = interp1d(distance_grid, integral_grid, kind="linear", bounds_error=False, fill_value=0.0)


#           matrix_element_fit = UnivariateSpline(distance_grid, integral_grid)

           # http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.interp1d.html
           # -> bounds_error and fill_value are flags we can control the value to be set if
           #    it is out of interpolation range. In this case, we set it to be 0.0.

           #print "test", 0.0, matrix_element_fit(0.5)
           #print "n_step"
           #print self.n_step
           #print len(integral_grid)

        return matrix_element_fit
    
    def slako_integrals_fit_from_C(self, distance_grid, dist_step, integral_grid, dist_atoms):

        n_grid_points = len(distance_grid)
        
        for i in xrange(n_grid_points):
            if integral_grid[i] != 0.0:
                inferior_integral_cutoff = distance_grid[i]
                index_inferior_cutoff = i
                break

        superior_integral_cutoff  = distance_grid[-1]
        index_superior_cutoff = n_grid_points

#         print "integral_cutoff ", inferior_integral_cutoff, superior_integral_cutoff
#         print "integral_cutoff ", index_inferior_cutoff, index_superior_cutoff

        initial_grid_point = 0
        final_grid_point = 0

        line_dist = int(round((dist_atoms/dist_step) - 1))

#         print n_grid_points, line_dist

        if dist_atoms >= inferior_integral_cutoff and dist_atoms <= (inferior_integral_cutoff + dist_step*2.0):
        
           initial_grid_point = index_inferior_cutoff
           final_grid_point = index_inferior_cutoff + 5

        elif dist_atoms >= (superior_integral_cutoff - dist_step*2.0) and dist_atoms <= superior_integral_cutoff:

           initial_grid_point = index_superior_cutoff - 5
           final_grid_point = index_superior_cutoff

        elif dist_atoms > superior_integral_cutoff:
           
	   matrix_element_fit = 0.0
	   return matrix_element_fit

        else:
            
           initial_grid_point = (line_dist - 2)
           final_grid_point = (line_dist + 3)

        # defining the ctypes:
        
        dimension = final_grid_point - initial_grid_point
#         print "dimension", dimension

        _distance_grid = c_double*(dimension)
        _distance_grid = _distance_grid()
        _integral_grid = c_double*(dimension)
        _integral_grid = _integral_grid()
        
        _dist_atoms = c_double(dist_atoms)
        _dist_atoms_fit = c_double(0.0)
        _err = c_double(0.0)
        
        count = 0
        for i in xrange(initial_grid_point, final_grid_point):

            _distance_grid[count] = distance_grid[i]
            _integral_grid[count] = integral_grid[i]

#            print distance_grid[i], integral_grid[i]
#            print _distance_grid[count], _integral_grid[count], count
            count = count + 1

        _order_of_fit = c_int(5)

        lib_polint.polint(_distance_grid,
                           _integral_grid,
                           _order_of_fit,
                           _dist_atoms,
                           byref(_dist_atoms_fit),
                           byref(_err))

        matrix_element_fit = float(_dist_atoms_fit.value)
 
        return matrix_element_fit

// This file is part of DFTBpy.
//
//    DFTBpy is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    DFTBpy is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with DFTBpy. If not, see <http://www.gnu.org/licenses/>.

// Created by MPL on 23rd of October 2014

# include <stdio.h>
# include <math.h>

void build_dgamma_matrix(int n_atoms, double dist_matrix[][n_atoms], char xyz_atom_symbols[n_atoms][3], double hubbard_s[n_atoms], double dgamma[][n_atoms]){

   int i, j, k;
   double ta, tb, fac, dfac, gabp, gbap, gab, gba, dgab, dgba;

   for (i = 0; i < n_atoms; i++){
       for (j = i; j < n_atoms; j++){

           ta = 3.2*hubbard_s[i];
           tb = 3.2*hubbard_s[j];

//           if (dist_matrix[i][j] == 0.0){
//               gamma[i][j] = 0.5*(( (ta*tb/(ta + tb))) + ( pow((ta*tb), 2)/pow((ta + tb), 3)) );
//           }

           if ( (dist_matrix[i][j] != 0.0) && (fabsf(ta - tb) < pow(10.0, -4)) ){
                 fac = ( (1.6*dist_matrix[i][j]*ta*tb)/(ta + tb) )*(1.0 + (ta*tb)/pow((ta + tb), 2) );
                 dfac = ( (1.6*ta*tb)/((ta + tb)*dist_matrix[i][j]) )*(1.0 + (ta*tb)/pow((ta + tb), 2) );
                 dgamma[i][j] = -1.0/pow(dist_matrix[i][j], 3) - (33.0*dfac + dfac*pow(fac, 2) + (9.0 + fac)*(2.0*fac*dfac))* \
                                (exp(-fac)/(48.0*dist_matrix[i][j])) - \
                                (48.0 + 33.0*fac + (9.0 + fac)*pow(fac, 2))*((-dfac*exp(-fac))/(48.0*dist_matrix[i][j]) - \
                                (exp(-fac)/(48.0*pow(dist_matrix[i][j], 3))));
           }

           if ( strcmp(xyz_atom_symbols[i], xyz_atom_symbols[j]) != 0 ) {

                dgamma[i][j] = -1.0/pow(dist_matrix[i][j], 3) + \
                               ( (-exp(-ta*dist_matrix[i][j]))* \
                               ((pow(tb, 6) -3.0*pow(tb, 4)*pow(ta, 2))/(pow(pow(ta, 2) - pow(tb, 2), 3)*pow(dist_matrix[i][j], 3))) - \
                               ((-ta*exp(-ta*dist_matrix[i][j]))/dist_matrix[i][j])* \
                               (((pow(tb, 4)*ta)/(2*pow( (pow(ta, 2) - pow(tb, 2)), 2))) - \
                               ((pow(tb, 6) - 3.0*pow(tb, 4)*pow(ta, 2))/(pow((pow(ta, 2) - pow(tb, 2)), 3)*dist_matrix[i][j])))) + \
                               ( (-exp(-tb*dist_matrix[i][j]))* \
                               ((pow(ta, 6) -3.0*pow(ta, 4)*pow(tb, 2))/(pow(pow(tb, 2) - pow(ta, 2), 3)*pow(dist_matrix[i][j], 3))) - \
                               ((-tb*exp(-tb*dist_matrix[i][j]))/dist_matrix[i][j])* \
                               (((pow(ta, 4)*tb)/(2*pow( (pow(tb, 2) - pow(ta, 2)), 2))) - \
                               ((pow(ta, 6) - 3.0*pow(ta, 4)*pow(tb, 2))/(pow((pow(tb, 2) - pow(ta, 2)), 3)*dist_matrix[i][j]))));



 
//                gabp = (pow(tb, 6) - 3.0*pow(ta, 2)*pow(tb, 4))/(dist_matrix[i][j]*pow( (pow(ta, 2) - pow(tb, 2) ), 3));
//                gbap = (pow(ta, 6) - 3.0*pow(tb, 2)*pow(ta, 4))/(dist_matrix[i][j]*pow( (pow(tb, 2) - pow(ta, 2) ), 3));
//                gab = exp(-ta*dist_matrix[i][j])*((0.5*ta*pow(tb, 4))/pow( (pow(ta, 2) - pow(tb, 2)), 2) - gabp);
//                gba = exp(-tb*dist_matrix[i][j])*((0.5*tb*pow(ta, 4))/pow( (pow(tb, 2) - pow(ta, 2)), 2) - gbap);
//                dgab = -ta*gab + (exp(-ta*dist_matrix[i][j])*gabp)/dist_matrix[i][j];
//                dgba = -tb*gba + (exp(-tb*dist_matrix[i][j])*gbap)/dist_matrix[i][j];
//                dgamma[i][j] = -(dgab+ dgba)/dist_matrix[i][j];
           }

           dgamma[j][i] = dgamma[i][j];
      }
   }

}

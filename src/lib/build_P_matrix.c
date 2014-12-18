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

// Created by MPL on 17th of October 2014

# include <stdio.h>

void build_p_matrix(int n_electrons, int n_total_orb, double eigen_vec[][n_total_orb], double eigen_val[n_total_orb], double P_dftb[][n_total_orb], double P_dftb_e[][n_total_orb]){

   int i, j, k;
   int n_closed_shell_electrons;
   n_closed_shell_electrons = n_electrons/2;

   for (i = 0; i < n_total_orb; i++){
       for (j = i; j < n_total_orb; j++){
            for (k = 0; k < n_closed_shell_electrons; k++){
                 P_dftb[i][j] = P_dftb[i][j] + 2.0*eigen_vec[i][k]*eigen_vec[j][k];
                 P_dftb_e[i][j] = P_dftb_e[i][j] + 2.0*eigen_val[k]*eigen_vec[i][k]*eigen_vec[j][k];
           }
           P_dftb[j][i] = P_dftb[i][j];
           P_dftb_e[j][i] = P_dftb_e[i][j];
       }
   }

}

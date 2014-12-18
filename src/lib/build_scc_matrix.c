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

void build_scc_matrix(int n_atoms, int n_total_orb, int atom_numb[n_total_orb], double S_dftb[][n_total_orb], double net_charge[n_atoms], double gamma[][n_atoms], double scc_matrix[][n_total_orb], double scc_matrix_without_S[][n_total_orb]){

   int i, j, k;

   for (i = 0; i < n_total_orb; i++){
       for (j = i; j < n_total_orb; j++){
            for (k = 0; k < n_atoms; k++){
                 scc_matrix[i][j] = scc_matrix[i][j] + 0.5*S_dftb[i][j]*(gamma[atom_numb[i]][k] + gamma[atom_numb[j]][k])*(-net_charge[k]);
                 scc_matrix_without_S[i][j] = scc_matrix_without_S[i][j] + 0.5*(gamma[atom_numb[i]][k] + gamma[atom_numb[j]][k])*(-net_charge[k]);
           }
           scc_matrix[j][i] = scc_matrix[i][j];
           scc_matrix_without_S[j][i] = scc_matrix_without_S[i][j];
       }
   }

}

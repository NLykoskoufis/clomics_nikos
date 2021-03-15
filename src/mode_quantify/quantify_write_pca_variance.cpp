/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


#include "quantify_data.h"

void quantify_data::writePhenotypesGroupsPCAVariance(string fvar)
{
    vrb.title("Writing proportion of variance in [" + fvar + "]");
    output_file fdo (fvar);
    
    fdo << "groups PC1 PC2" << std::endl;
    for(int i = 0; i < vec_groups.size(); i++)
    {
        fdo << phenotype_id[vec_groups[i][0]] << " " << phenotype_id[vec_groups[i][1]];
        for(int p =0; p < vec_prop_var[i].size(); p++)
        {
            fdo << " " << vec_prop_var[i][p];
        }
        fdo << std::endl;
    }
}

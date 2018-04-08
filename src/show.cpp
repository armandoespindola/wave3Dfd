/*
                  ### Wave3Dfd ####

    Copyright (C) April 2018  Armando Espindola Carmona,
    Universidad Nacional Autonoma de Mexico (UNAM)
    King Abdullah University of Science and Technology (KAUST).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "show.hpp"

void Show::print(Dfloat f){
	 printf("%f\n",f);
}

void Show::print(Dfloat *f, int l){

	for (int i = 0; i<l; ++i)
	{

	  printf("%f\t%d\n",f[i],i);

	}

}

void Show::print(Dfloat *f, int l, int m){

	for (int j = 0; j<m; ++j){
		for (int i = 0; i<l; ++i){


		 printf("%f\t%d\t%d\n",f[i + j * l],i,j);



		}
	}

}
	

		
void Show::print(Dfloat *f, int l, int m, int n){

for (int k = 0; k<n; ++k){
	for (int j = 0; j<m; ++j){
		for (int i = 0; i<l; ++i){


			printf("  %f",f[i + j * l + k * l * m]);



		}

		printf("\n");


	        
	}

	printf("\n");
	printf("\n");

}


}	

		
	
	

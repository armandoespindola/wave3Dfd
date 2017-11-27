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


			printf("  %10.4f",f[i + j * l + k * l * m]);



		}

		printf("\n");


	        
	}

	printf("\n");
	printf("\n");

}


}	

		
	
	

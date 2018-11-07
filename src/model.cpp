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

#include "model.hpp"

MODEL::MODEL(std::string FileP,std::string FileS,std::string FileR, VecI iGDim, VecI iSubDomNodeN) {

  GDim = iGDim;
  SubDomNodeN = iSubDomNodeN;
  FILE *R,*P,*S;
  int nsort;

  R=fopen(FileR.c_str(),"rb");
  P=fopen(FileP.c_str(),"rb");
  S=fopen(FileS.c_str(),"rb");

  NDT.x = GDim.x + 2 * PML.x;
  NDT.y = GDim.y + 2 * PML.y;
  NDT.z = GDim.z +  PML.z;

  NDL.x = SubDomNodeN.x + 2 * KHALO;
  NDL.y = SubDomNodeN.y + 2 * KHALO;
  NDL.z = SubDomNodeN.z + 2 * KHALO;

  mu = new Dfloat[NDT.x * NDT.y * NDT.z];
  lamb = new Dfloat[NDT.x * NDT.y * NDT.z];
  rho = new Dfloat[NDT.x * NDT.y * NDT.z];

  muL = new Dfloat[NDL.x * NDL.y * NDL.z];
  lambL = new Dfloat[NDL.x * NDL.y * NDL.z];
  rhoL = new Dfloat[NDL.x * NDL.y * NDL.z];


  Nsub.x = NDT.x / SubDomNodeN.x;
  Nsub.y = NDT.y / SubDomNodeN.y;
  Nsub.z = NDT.z / SubDomNodeN.z;  

  Dfloat rho_buff,vp_buff,vs_buff;
  float *Rbuff,*Sbuff,*Pbuff;

  Rbuff = new float[GDim.x*GDim.y*GDim.z];
  Sbuff = new float[GDim.x*GDim.y*GDim.z];
  Pbuff = new float[GDim.x*GDim.y*GDim.z];

  fread(Rbuff,sizeof(float),GDim.x*GDim.y*GDim.z*sizeof(float),R);
  fread(Pbuff,sizeof(float),GDim.x*GDim.y*GDim.z*sizeof(float),P);
  fread(Sbuff,sizeof(float),GDim.x*GDim.y*GDim.z*sizeof(float),S);

  for (int k=0;k<GDim.z;k++){
    for (int j=0;j<GDim.y;j++){
      for (int i=0;i<GDim.x;i++){

        int indx =  (i + PML.x) + (j + PML.y) * NDT.x + (k+ PML.z) * NDT.x * NDT.y;
        int idxL = i + j * GDim.x + k * GDim.x * GDim.y;
        rho[indx] = (Dfloat) Rbuff[idxL];
        mu[indx] = pow((Dfloat) Sbuff[idxL],2.0) * rho[indx];
        lamb[indx]= rho[indx]* pow((Dfloat) Pbuff[idxL],2.0) - (2.0 * mu[indx]);

      }
    }
  }


 delete [] Rbuff;
 delete [] Sbuff;
 delete [] Pbuff;


  fclose(R);
  fclose(P);
  fclose(S);

  PML_MODEL();


}


MODEL::~MODEL(){

  delete [] mu;
  delete [] rho;
  delete [] lamb;
  delete [] muL;
  delete [] rhoL;
  delete [] lambL;

}




void MODEL::PML_MODEL(){

  int indx1,indx2;

  // X DIRECTION
  
  for (int k=0;k<GDim.z;k++){
    for (int j=0;j<GDim.y;j++){
      for (int i=0;i<PML.x;i++){

	indx1 = i + (PML.y + j) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	indx2 = PML.x + (j + PML.y) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;

	// LEFT SIDE
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];

	// RIGHT SIDE
	

	indx1 = (PML.x + GDim.x + i) + (j+PML.y) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	
  indx2 = (PML.x + GDim.x - 1) + (j+PML.y) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];
				
      }
    }
  }


  // Y DIRECTION
  
    for (int k=0;k<GDim.z;k++){
    for (int j=0;j<PML.y;j++){
      for (int i=0;i<NDT.x;i++){

	indx1 = i + j * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	indx2 = i + PML.y * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;

	// LEFT SIDE
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];

	// RIGHT SIDE
	

	indx1 = i + (PML.y + GDim.y + j) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	indx2 = i + (PML.y + GDim.y - 1) * NDT.x + \
  (PML.z + k) * NDT.x * NDT.y;
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];
				
      }
    }
  }

      // Z DIRECTION
  
    for (int k=0;k<PML.z;k++){
    for (int j=0;j<NDT.y;j++){
      for (int i=0;i<NDT.x;i++){

	indx1 = i + j * NDT.x + k * NDT.x * NDT.y;
	indx2 = i + j * NDT.x + PML.z * NDT.x * NDT.y;

	// LEFT SIDE
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];
				
      }
    }
  }


}

void MODEL::SortModel(VecI subi[], Dfloat *SortMod){

  int indx1,indx2;
  int ii,jj,kk;
  // Local Subdomain Nodes
  int nodl = NDL.x * NDL.y * NDL.z; 


 for (int idx=0;idx<Nsub.x * Nsub.y * Nsub.z;idx++){
	 //std::cout<<subi[idx].x<<subi[idx].y<<subi[idx].z<<std::endl;
  for (int k=0;k<NDL.z;k++){
    for (int j=0;j<NDL.y;j++){
      for (int i=0;i<NDL.x;i++){

	ii = i + (subi[idx].x * SubDomNodeN.x) - KHALO;
	jj = j + (subi[idx].y * SubDomNodeN.y) - KHALO;
	kk = k + (subi[idx].z * SubDomNodeN.z) - KHALO;


	indx1 = i + j * NDL.x + k * NDL.x * NDL.y;
	indx2 = ii + jj * NDT.x + kk * NDT.x * NDT.y;

	if ( (ii < 0) || (jj < 0) || (kk < 0) || \
	     (ii >= NDT.x) || (jj >= NDT.y) || (kk >= NDT.z) ) {
	  
	//SubR[indx1] = 0.0;
	SortMod[indx1 + idx * 3 * nodl] = 0.0;
	//SubM[indx1] = 0.0;
	SortMod[indx1 + nodl + idx * 3 * nodl] = 0.0;
	//SubL[indx1] = 0.0;
	SortMod[indx1 + 2 * nodl + idx * 3 * nodl] = 0.0;

	} else {

	  //SubR[indx1] = rho[indx2];
	  SortMod[indx1 + idx * 3 * nodl] = rho[indx2];
	  //SubM[indx1] = mu[indx2];
	  SortMod[indx1 + nodl + idx * 3 * nodl] = mu[indx2];
	  //SubL[indx1] = lamb[indx2];
	  SortMod[indx1 + 2 * nodl + idx * 3 * nodl] = lamb[indx2];
	
	}
	

      }
    }
  }

 }
  
  
}










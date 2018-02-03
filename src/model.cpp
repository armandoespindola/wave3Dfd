#include "model.hpp"

MODEL::MODEL(std::string FileP,std::string FileS,std::string FileR, VecI iGDim, VecI iSubDomNodeN) {

  GDim = iGDim;
  SubDomNodeN = iSubDomNodeN;
  FILE *R,*P,*S;

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
  float buff;

  for (int k=0;k<GDim.z;k++){
    for (int j=0;j<GDim.y;j++){
      for (int i=0;i<GDim.x;i++){

	int indx =  (i + PML.x) + (j + PML.y) * NDT.x + (k+ PML.z) * NDT.x * NDT.y; 
  
	fread(&buff,sizeof(float),1,R);
  rho[indx] = (Dfloat) buff;
	fread(&buff,sizeof(float),1,P);
  vp_buff = (Dfloat) buff;
	fread(&buff,sizeof(float),1,S);
  vs_buff = (Dfloat) buff;


	mu[indx] = pow(vs_buff,2.0) * rho[indx];
	lamb[indx]= rho[indx]* pow(vp_buff,2.0) - (2.0 * mu[indx]);
	
      }
    }
  }


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

void MODEL::SubModel(VecI NumSub, Dfloat *SubR, Dfloat *SubM, Dfloat *SubL){

  int indx1,indx2;
  int ii,jj,kk;

  for (int k=0;k<NDL.z;k++){
    for (int j=0;j<NDL.y;j++){
      for (int i=0;i<NDL.x;i++){

	ii = i + (NumSub.x * SubDomNodeN.x) - KHALO;
	jj = j + (NumSub.y * SubDomNodeN.y) - KHALO;
	kk = k + (NumSub.z * SubDomNodeN.z) - KHALO;


	indx1 = i + j * NDL.x + k * NDL.x * NDL.y;
	indx2 = ii + jj * NDT.x + kk * NDT.x * NDT.y;

	if ( (ii < 0) || (jj < 0) || (kk < 0) || \
	     (ii >= NDT.x) || (jj >= NDT.y) || (kk >= NDT.z) ) {
	  
	SubR[indx1] = 0.0;
	SubM[indx1] = 0.0;
	SubL[indx1] = 0.0;

	} else {

	  SubR[indx1] = rho[indx2];
	  SubM[indx1] = mu[indx2];
	  SubL[indx1] = lamb[indx2];
	
	}
	

      }
    }
  }

  
  
}










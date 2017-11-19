#include "model.hpp"

MODEL::MODEL(char *FileP,char *FileS,char *FileR, VecI iGDim, VecI iSubDomNodeN) {

  GDim = iGDim;
  SubDomNodeN = iSubDomNodeN;
  FILE *R,*P,*S;

  R=fopen(FileR,"rb");
  P=fopen(FileP,"rb");
  S=fopen(FileS,"rb");

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

  for (int k=0;k<GDim.z;k++){
    for (int j=0;j<GDim.y;j++){
      for (int i=0;i<GDim.x;i++){

	int indx =  (i + PML.x) + (j + PML.y) * NDT.x + (k+ PML.z) * NDT.x * NDT.y; 
  
	fread(&rho[indx],sizeof(Dfloat),1,R);
	fread(&vp_buff,sizeof(Dfloat),1,P);
	fread(&vs_buff,sizeof(Dfloat),1,S);

	mu[indx] = pow(vs_buff,2.0) * rho[indx];
	lamb[indx]= rho[indx]* pow(vp_buff,2.0) - (2.0 * mu[indx]);
	
      }
    }
  }


  fclose(R);
  fclose(P);
  fclose(S);

  //PML_MODEL();
  
  


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

	indx1 = i + j * NDT.x + k * NDT.x * NDT.y;
	indx2 = PML.x + j * NDT.x + k * NDT.x * NDT.y;

	// LEFT SIDE
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];

	// RIGHT SIDE
	

	indx1 = (PML.x + NDT.x + i) + j * NDT.x + k * NDT.x * NDT.y;
	indx2 = (PML.x + NDT.x - 1) + j * NDT.x + k * NDT.x * NDT.y;
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];
				
      }
    }
  }


  // Y DIRECTION
  
    for (int k=0;k<GDim.z;k++){
    for (int j=0;j<PML.y;j++){
      for (int i=0;i<GDim.x;i++){

	indx1 = i + j * NDT.x + k * NDT.x * NDT.y;
	indx2 = i + PML.y * NDT.x + k * NDT.x * NDT.y;

	// LEFT SIDE
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];

	// RIGHT SIDE
	

	indx1 = i + (PML.x + NDT.x + j) * NDT.x + k * NDT.x * NDT.y;
	indx2 = i + (PML.x + NDT.x - 1) * NDT.x + k * NDT.x * NDT.y;
	
	mu[indx1] = mu[indx2];
	rho[indx1] = rho[indx2];
	lamb[indx1] = lamb[indx2];
				
      }
    }
  }

      // Z DIRECTION
  
    for (int k=0;k<PML.z;k++){
    for (int j=0;j<GDim.y;j++){
      for (int i=0;i<GDim.x;i++){

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

	ii = i + NumSub.x * SubDomNodeN.x;
	jj = j + NumSub.y * SubDomNodeN.y;
	kk = k + NumSub.z * SubDomNodeN.z;


	indx1 = i + j * NDL.y + k * NDL.x * NDL.y;
	
	indx2 = ii + jj * NDT.y + kk * NDT.x * NDT.y;

	if ( ((ii-KHALO) < 0) || ((jj-KHALO) < 0) || ((kk-KHALO) < 0) || \
	     ((ii+KHALO) >= NDT.x) || ((jj+KHALO) >= NDT.y) || ((kk+KHALO) >= NDT.z) ) {
	  
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










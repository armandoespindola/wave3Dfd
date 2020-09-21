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

#include "sdm.hpp"

SDM::SDM(VecF IGI, VecF IGF,VecI I_NodG,VecF IlimI, VecF IlimF, VecI INod, \
	 VecI ICPML,Dfloat If0, Dfloat Idt, VecI INsdm, VecI INumSubDom,int iPROPAGATION) {

	// GI INITIAL GLOBAL LIMIT WITH PML 
	// GF END GLOBAL LIMIT WITH PML
	// NodG NUMBER OF GLOBAL NODES WITHOUT PML
	// IlimI INITIAL LIMIT SUB DOMAIN 
	// IlimF END LIMIT SUB DOMAIN
	// INOD   NUMBER OF LOCAL NODES WITHOUT HALO SUBDOMAIN
        // ICPML PML NODES
	// f0 FREQUENCY
	// dt DELTA T
        // Nsdm SUBDOMAIN INDEX
  // NumSubDom Number total of subdomains

	GI = IGI;
	GF = IGF;
	f0 = If0;
	dt = Idt;
	PML = ICPML;
	NodG = I_NodG;
	HALO.x = KHALO;
	HALO.y = KHALO;
	HALO.z = KHALO;
	NodLoc = INod;
	Nsdm = INsdm;
	NumSubDom = INumSubDom;
	PROPAGATION = iPROPAGATION;

	if (PROPAGATION == 0){
	  sgn = 1.0;
	}else if (PROPAGATION == 1 ){
	  sgn = -1.0;
	}
	
	VecI ELE={INod.x-1,INod.y-1,INod.z-1};

	SDMGeom = new geometry3D(IlimI,IlimF,ELE,HALO);

	thickness_PML = {SDMGeom->Dx() * PML.x,SDMGeom->Dy() * PML.y,SDMGeom->Dz() * PML.z};

	NT.x = NodG.x + 2 * PML.x;
	NT.y = NodG.y + 2 * PML.y;
	NT.z = NodG.z + PML.z;

	pml_x = new DPML(NodG.x,NT.x,thickness_PML.x,SDMGeom->Dx(),PML.x,dt,f0, \
		PML_XMIN,PML_XMAX);

	pml_y = new DPML(NodG.y,NT.y,thickness_PML.y,SDMGeom->Dy(),PML.y,dt,f0, \
		PML_YMIN,PML_YMAX);

	pml_z = new DPML(NodG.z,NT.z,thickness_PML.z,SDMGeom->Dz(),PML.z,dt,f0, \
		PML_ZMIN,PML_ZMAX);



	sxx = new Dfloat [SDMGeom->HALO_Node()];
	syy = new Dfloat [SDMGeom->HALO_Node()];
	szz = new Dfloat [SDMGeom->HALO_Node()];
	sxy = new Dfloat [SDMGeom->HALO_Node()];
	sxz = new Dfloat [SDMGeom->HALO_Node()];
	syz = new Dfloat [SDMGeom->HALO_Node()];
	vx = new Dfloat [SDMGeom->HALO_Node()];
	vy = new Dfloat [SDMGeom->HALO_Node()];
	vz = new Dfloat [SDMGeom->HALO_Node()];
	ux = new Dfloat [SDMGeom->HALO_Node()];
	uy = new Dfloat [SDMGeom->HALO_Node()];
	uz = new Dfloat [SDMGeom->HALO_Node()];
	mu = new Dfloat [SDMGeom->HALO_Node()];
	lamb = new Dfloat [SDMGeom->HALO_Node()];
	rho = new Dfloat [SDMGeom->HALO_Node()];
	dsxx_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsxy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dsxz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dsxy_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsyy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dsyz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dsxz_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsyz_dy = new Dfloat [SDMGeom->HALO_Node()];
	dszz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dy = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dy = new Dfloat [SDMGeom->HALO_Node()];







	if (Nsdm.x == 0){
	  // LEFT X
	  bn_lx = new Dfloat[SDMGeom->L_NodeZ() * SDMGeom->L_NodeY() * KHALO * 9 * nsteps];
	}else if (Nsdm.x == NumSubDom.x-1) {
	  // RUGHT X
	  bn_rx = new Dfloat[SDMGeom->L_NodeZ() * SDMGeom->L_NodeY() * KHALO * 9 * nsteps];
	}
	if (Nsdm.y == 0){
	  // LEFT Y
	  bn_ly = new Dfloat[SDMGeom->L_NodeZ() * SDMGeom->L_NodeX() * KHALO * 9 * nsteps];
	}else if (Nsdm.y == NumSubDom.y-1) {
	  // RIGHT Y
	  bn_ry = new Dfloat[SDMGeom->L_NodeZ() * SDMGeom->L_NodeX() * KHALO * 9 * nsteps];
	}
	if (Nsdm.z == 0){
	  // LEFT Z
	  bn_lz = new Dfloat[SDMGeom->L_NodeY() * SDMGeom->L_NodeX() * KHALO * 9 * nsteps];
	}

}


SDM::~SDM(){

	delete  pml_x;
	delete  pml_y;
	delete  pml_z;
	delete  SDMGeom;
	delete [] sxx;
	delete [] syy;
	delete [] szz;
	delete [] sxy;
	delete [] sxz;
	delete [] syz;
	delete [] vx;
	delete [] vy;
	delete [] vz;
	delete [] ux;
	delete [] uy;
	delete [] uz;
	delete [] mu;
	delete [] lamb;
	delete [] rho;
	delete [] dsxx_dx;
	delete [] dsxy_dy;
	delete [] dsxz_dz;
	delete [] dsxy_dx;
	delete [] dsyy_dy;
	delete [] dsyz_dz;
	delete [] dsxz_dx;
	delete [] dsyz_dy;
	delete [] dszz_dz;
	delete [] dvx_dx;
	delete [] dvy_dy;
	delete [] dvz_dz;
	delete [] dvx_dy;
	delete [] dvy_dx;
	delete [] dvx_dz;
	delete [] dvz_dx;
	delete [] dvy_dz;
	delete [] dvz_dy;


	if (Nsdm.x == 0){
	  // LEFT X
	  delete [] bn_lx;
	}else if (Nsdm.x == NumSubDom.x-1) {
	  // RUGHT X
	  delete [] bn_rx;
	}
	if (Nsdm.y == 0){
	  // LEFT Y
	  delete [] bn_ly;
	}else if (Nsdm.y == NumSubDom.y-1) {
	  // RIGHT Y
	  delete [] bn_ry;
	}
	if (Nsdm.z == 0){
	  // LEFT Z
	  delete [] bn_lz;
	}

}


// INDEX_IJK

int SDM::IJK(int i, int j, int k){

  int indx = i + j * SDMGeom->HALO_NodeX() + k * SDMGeom->HALO_NodeX() * \
    SDMGeom->HALO_NodeY();

  return indx;

}


void SDM::ModelRead(Dfloat *model,char param[]){

  if (strcmp("RHO",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i){
      rho[i] = model[i];
    }


  }


  if (strcmp("MU",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i){
      mu[i] = model[i];


    }

  }


  if (strcmp("LAMB",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i){
      lamb[i] = model[i];
    }

  }

}



void SDM::InitVar(Dfloat f){


  for (int i=0; i<SDMGeom->HALO_Node(); ++i){


    sxx[i] = f;
    syy[i] = f;
    szz[i] = f;
    sxy[i] = f;
    sxz[i] = f;
    syz[i] = f;
    vx[i] = f;
    vy[i] = f;
    vz[i] = f;
    ux[i] = f;
    uy[i] = f;
    uz[i] = f;
    dsxx_dx[i] = f;
    dsxy_dy[i] = f;
    dsxz_dz[i] = f;
    dsxy_dx[i] = f;
    dsyy_dy[i] = f;
    dsyz_dz[i] = f;
    dsxz_dx[i] = f;
    dsyz_dy[i] = f;
    dszz_dz[i] = f;
    dvx_dx[i] = f;
    dvy_dy[i] = f;
    dvz_dz[i] = f;
    dvx_dy[i] = f;
    dvy_dx[i] = f;
    dvx_dz[i] = f;
    dvz_dx[i] = f;
    dvy_dz[i] = f;
    dvz_dy[i] = f;

  }

  

}


void SDM::EB(Dfloat *BN,Dfloat* DomLoc, char *param){
  int indx1,indx2;
  
  if (strcmp("South",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<HALO.y;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * HALO.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }


  }else if (strcmp("North",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<HALO.y;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * HALO.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + SDMGeom->HALO_NodeY() - 2*HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }


  }else if (strcmp("West",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<HALO.x;i++){

	  indx1 = i + j * HALO.x + k * NodLoc.y * HALO.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }




  }else if (strcmp("East",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<HALO.x;i++){

	  indx1 = i + j * HALO.x + k * NodLoc.y * HALO.x;
	  indx2 = (SDMGeom->HALO_NodeX() + i - 2*HALO.x) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }



  }else if (strcmp("UP",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<HALO.z;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * NodLoc.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + SDMGeom->HALO_NodeZ() - 2*HALO.z ) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }



  }else if (strcmp("DOWN",param) == 0){


#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<HALO.z;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * NodLoc.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z ) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();	  

	  BN[indx1] = DomLoc[indx2];
	  
	}
      }
    }



  } else {

    std::cout<<"SDM::EB THAT IS NOT A BOUNDARY OPTION:"<<param<<std::endl;


  }


}


void SDM::IB(Dfloat *BN, Dfloat *DomLoc, char *param) {
  int indx1,indx2;

  if (strcmp("South",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<HALO.y;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * HALO.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }


  }else if (strcmp("North",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<HALO.y;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * HALO.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + SDMGeom->HALO_NodeY() - HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }


  }else if (strcmp("West",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<HALO.x;i++){

	  indx1 = i + j * HALO.x + k * NodLoc.y * HALO.x;
	  indx2 = (i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }




  }else if (strcmp("East",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<HALO.x;i++){

	  indx1 = i + j * HALO.x + k * NodLoc.y * HALO.x;
	  indx2 = (SDMGeom->HALO_NodeX() + i - HALO.x) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY() ;	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }



  }else if (strcmp("UP",param) == 0){

#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<HALO.z;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * NodLoc.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k + SDMGeom->HALO_NodeZ() - HALO.z ) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }



  }else if (strcmp("DOWN",param) == 0){


#pragma omp parallel for num_threads(N_omp)\
  private(indx1,indx2)
    for (int k=0;k<HALO.z;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx1 = i + j * NodLoc.x + k * NodLoc.y * NodLoc.x;
	  indx2 = (HALO.x + i) + (j + HALO.y) * SDMGeom->HALO_NodeX() \
	    + (k) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();	  

	  DomLoc[indx2] = BN[indx1];
	  
	}
      }
    }



  } else {

    std::cout<<"SDM::IB THAT IS NOT A BOUNDARY OPTION:"<<param<<std::endl;


  }


}



void SDM::ExpBoundary(Dfloat *BN, char *TBound,char *NameVar) {

  if (strcmp("SXX",NameVar) == 0){
    EB(BN,sxx,TBound);
  } else if (strcmp("SYY",NameVar) == 0){
    EB(BN,syy,TBound);
  }else if (strcmp("SZZ",NameVar) == 0){
    EB(BN,szz,TBound);
  } else if (strcmp("SXY",NameVar) == 0){
    EB(BN,sxy,TBound);
  }else if (strcmp("SXZ",NameVar) == 0){
    EB(BN,sxz,TBound);
  } else if (strcmp("SYZ",NameVar) == 0){
    EB(BN,syz,TBound);
  }else if (strcmp("VX",NameVar) == 0){
    EB(BN,vx,TBound);
  }else if (strcmp("VY",NameVar) == 0){
    EB(BN,vy,TBound);
  }else if (strcmp("VZ",NameVar) == 0){
    EB(BN,vz,TBound);
  }else if (strcmp("UX",NameVar) == 0){
    EB(BN,ux,TBound);
  }else if (strcmp("UY",NameVar) == 0){
    EB(BN,uy,TBound);
  }else if (strcmp("UZ",NameVar) == 0){
    EB(BN,uz,TBound);
  } else {
    std::cout<<"SDM::ExpBoundary:: IT IS NOT A BOUNDARY OPTION:"<<NameVar<<std::endl;
  }
  
}


void SDM::ImpBoundary(Dfloat *BN, char *TBound,char *NameVar) {

  if (strcmp("SXX",NameVar) == 0){
    IB(BN,sxx,TBound);
  } else if (strcmp("SYY",NameVar) == 0){
    IB(BN,syy,TBound);
  } else if (strcmp("SZZ",NameVar) == 0){
    IB(BN,szz,TBound);
  } else if (strcmp("SXY",NameVar) == 0){
    IB(BN,sxy,TBound);
  } else if (strcmp("SXZ",NameVar) == 0){
    IB(BN,sxz,TBound);
  } else if (strcmp("SYZ",NameVar) == 0){
    IB(BN,syz,TBound);
  } else if (strcmp("VX",NameVar) == 0){
    IB(BN,vx,TBound);
  } else if (strcmp("VY",NameVar) == 0){
    IB(BN,vy,TBound);
  }else if (strcmp("VZ",NameVar) == 0){
    IB(BN,vz,TBound);
  } else if (strcmp("UX",NameVar) == 0){
    IB(BN,ux,TBound);
  } else if (strcmp("UY",NameVar) == 0){
    IB(BN,uy,TBound);
  }else if (strcmp("UZ",NameVar) == 0){
    IB(BN,uz,TBound);
  } else {
    std::cout<<"SDM::ImpBoundary:: IT IS NOT A BOUNDARY OPTION:"<<NameVar<<std::endl;
  }
}

Dfloat SDM::SCoorX(int indx){
  Dfloat coor = SDMGeom->CoorX(indx) - SDMGeom->thickness_PML().x;
  return coor;
}

Dfloat SDM::SCoorY(int indx){
  Dfloat coor = SDMGeom->CoorY(indx) - SDMGeom->thickness_PML().y;
  return coor;
}

Dfloat SDM::SCoorZ(int indx){
  Dfloat coor = SDMGeom->CoorZ(indx) - SDMGeom->thickness_PML().z;
  return coor;
}

Dfloat SDM::SCoorXHalf(int indx){
  Dfloat coor = SDMGeom->CoorXHalf(indx) - SDMGeom->thickness_PML().x;
  return coor;
}

Dfloat SDM::SCoorYHalf(int indx){
  Dfloat coor = SDMGeom->CoorYHalf(indx) - SDMGeom->thickness_PML().y;
  return coor;
}

Dfloat SDM::SCoorZHalf(int indx){
  Dfloat coor = SDMGeom->CoorZHalf(indx) - SDMGeom->thickness_PML().z;
  return coor;
}


VecI SDM::Loc2Glo(VecI indx){
  VecI ind;
  
  ind.x = indx.x - HALO.x + Nsdm.x * NodLoc.x;
  ind.y = indx.y - HALO.y + Nsdm.y * NodLoc.y;
  ind.z = indx.z - HALO.z + Nsdm.z * NodLoc.z;
  
  return ind;

}


VecI SDM::SFindNode(VecI Nodef){

  VecI ind = {-1,-1,-1};

  VecI Global,Local;
  
  for (int iz=0;iz<SDMGeom->L_NodeZ();iz++){
    Local = {HALO.x,HALO.y,iz+HALO.z};
    Global = Loc2Glo(Local);

    if ((Nodef.z - Global.z) == 0){
      ind.z = iz;
    }

  }

  
  for (int iy=0;iy<SDMGeom->L_NodeY();iy++){
    Local = {HALO.x,iy+HALO.y,HALO.z};
    Global = Loc2Glo(Local);

    if ((Nodef.y - Global.y) == 0){
      ind.y = iy;
    }

  }


  for (int ix=0;ix<SDMGeom->L_NodeX();ix++){
    Local = {ix+HALO.x,HALO.y,HALO.z};
    Global = Loc2Glo(Local);

    if ((Nodef.x - Global.x) == 0){
      ind.x = ix;
    }

  }

 return ind;
}

void SDM::AddVal(VecI indx, char *NameVar, Dfloat Val){
  VecI ind;
  int i;
  
  ind = SFindNode(indx);

  if ((ind.x > -1) && (ind.y > -1) && (ind.z > -1)) {

   ind.x += HALO.x;
   ind.y += HALO.y;
   ind.z += HALO.z;

   if (strcmp("SXX",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    sxx[i] += Val;
  } else if (strcmp("SYY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    syy[i] += Val;
  } else if (strcmp("SZZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    szz[i] += Val;
  } else if (strcmp("SXY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    sxy[i] += Val;
  } else if (strcmp("SXZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    sxz[i] += Val;
  } else if (strcmp("SYZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    syz[i] += Val;
  } else if (strcmp("VX",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    vx[i] += Val;
  } else if (strcmp("VY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    vy[i] += Val;
  } else if (strcmp("VZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    vz[i] += Val;
  } else {
    std::cout<<"SDM::AddVal:: IT IS NOT A VARIABLE OPTION:"<<NameVar<<std::endl;
  }
}
  
}


Dfloat SDM::GetVal(VecI indx, char *NameVar){
  VecI ind;
  Dfloat Val;
  int i;

  ind = SFindNode(indx);

  if ((ind.x > -1) && (ind.y > -1) && (ind.z > -1)) {

   ind.x += HALO.x;
   ind.y += HALO.y;
   ind.z += HALO.z;
  
   if (strcmp("SXX",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = sxx[i];
  } else if (strcmp("SYY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = syy[i];
  } else if (strcmp("SZZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = szz[i];
  } else if (strcmp("SXY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = sxy[i];
  } else if (strcmp("SXZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = sxz[i];
  } else if (strcmp("SYZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = syz[i];
  } else if (strcmp("VX",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = vx[i];
  } else if (strcmp("VY",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = vy[i];
  } else if (strcmp("VZ",NameVar) == 0){
    i = ind.x + ind.y * SDMGeom->HALO_NodeX() + ind.z * SDMGeom->HALO_NodeX() * \
      SDMGeom->HALO_NodeY();
    Val = vz[i];
  } else {
    std::cout<<"SDM::GetVal:: IT IS NOT A VARIABLE OPTION:"<<NameVar<<std::endl;
  }
  return Val;
}

}

void SDM:: FD_SII(VecI Init,VecI Iend){
  VecI Lindx, Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg,lamb_avg;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,lamb_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (vx[IJK(ix+1,iy,iz)] - vx[IJK(ix,iy,iz)]) - \
		  C0 * (vx[IJK(ix+2,iy,iz)] - vx[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (vy[IJK(ix,iy,iz)] - vy[IJK(ix,iy-1,iz)]) - \
		  C0 * (vy[IJK(ix,iy+1,iz)] - vy[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();
	
	df_dK = ( C1 * (vz[IJK(ix,iy,iz)] - vz[IJK(ix,iy,iz-1)]) - \
		  C0 * (vz[IJK(ix,iy,iz+1)] - vz[IJK(ix,iy,iz-2)]) ) / SDMGeom->Dz();


	if (PROPAGATION == 0) {
	  
	dvx_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dvx_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dvy_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dvy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dvz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dvz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dvx_dx[IJK(ix,iy,iz)];
	df_dJ += dvy_dy[IJK(ix,iy,iz)];
	df_dK += dvz_dz[IJK(ix,iy,iz)];

	
	}
	
	// Average properties
	
	mu_avg =   1.0 / ((0.50 / mu[IJK(ix,iy,iz)]) + (0.50 / mu[IJK(ix+1,iy,iz)]));
	lamb_avg = 1.0 / ((0.50 / lamb[IJK(ix,iy,iz)]) + (0.50 / lamb[IJK(ix+1,iy,iz)]));
	
	// SXX

	sxx[IJK(ix,iy,iz)] = sxx[IJK(ix,iy,iz)] + sgn * \
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dI) + sgn * \
	  dt * lamb_avg *  (df_dJ + df_dK);

	// SYY
	
	syy[IJK(ix,iy,iz)] = syy[IJK(ix,iy,iz)] + sgn *\
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dJ) + sgn *\
	  dt * lamb_avg *  (df_dI + df_dK);

	// SZZ
	
	szz[IJK(ix,iy,iz)] = szz[IJK(ix,iy,iz)] + sgn * \
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dK) +  sgn *\
	  dt * lamb_avg *  (df_dJ + df_dI);


	

      }
    }
  }

}


void SDM::Free_VX(VecI Init,VecI Iend,int zh){
  int iz = Iend.z + HALO.z - zh ;
  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,Lindx,Gindx)\
  firstprivate(zh,iz)
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxx[IJK(ix,iy,iz)] - sxx[IJK(ix-1,iy,iz)]) - \
		  C0 * (sxx[IJK(ix+1,iy,iz)] - sxx[IJK(ix-2,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (sxy[IJK(ix,iy,iz)] - sxy[IJK(ix,iy-1,iz)]) - \
		  C0 * (sxy[IJK(ix,iy+1,iz)] - sxy[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();


	if (zh == 1) { // Free Surface Z = 0;
	
	  df_dK = -( (35.0 / 8.0) * sxz[IJK(ix,iy,iz-1)] \
		    - (35.0 / 24.0) * sxz[IJK(ix,iy,iz-2)] \
		    + (21.0 / 40.0) * sxz[IJK(ix,iy,iz-3)] \
		    - (5.0 / 56.0) * sxz[IJK(ix,iy,iz-4)] ) / SDMGeom->Dz();

	} else if (zh == 2) { // Free Surface Z == h;


	  df_dK = -( - (31.0 / 24.0) * sxz[IJK(ix,iy,iz)] \
		    + (29.0 / 24.0) * sxz[IJK(ix,iy,iz-1)] \
		    - (3.0 / 40.0) * sxz[IJK(ix,iy,iz-2)] \
		    + (1.0 / 168.0) * sxz[IJK(ix,iy,iz-3)] ) / SDMGeom->Dz();
	  
	}


	if (PROPAGATION == 0) {
	  
	dsxx_dx[IJK(ix,iy,iz)] = pml_x->B_HALF(Gindx.x) * dsxx_dx[IJK(ix,iy,iz)] \
	  + pml_x->A_HALF(Gindx.x) * df_dI;

	dsxy_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dsxy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dsxz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dsxz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dsxx_dx[IJK(ix,iy,iz)];
	df_dJ += dsxy_dy[IJK(ix,iy,iz)];
	df_dK += dsxz_dz[IJK(ix,iy,iz)];

	}


	if (PROPAGATION == 1){
	  ux[IJK(ix,iy,iz)] = ux[IJK(ix,iy,iz)] + sgn * vx[IJK(ix,iy,iz)] * dt;
	}
	
	vx[IJK(ix,iy,iz)] = vx[IJK(ix,iy,iz)] + sgn * (dt / rho[IJK(ix,iy,iz)]) * \
	  (df_dI + df_dJ + df_dK);

	if (PROPAGATION == 0){
	  ux[IJK(ix,iy,iz)] = ux[IJK(ix,iy,iz)] + sgn * vx[IJK(ix,iy,iz)] * dt;
	  
	}

	
      }
    }

  
}


void SDM::FD_VX(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxx[IJK(ix,iy,iz)] - sxx[IJK(ix-1,iy,iz)]) - \
		  C0 * (sxx[IJK(ix+1,iy,iz)] - sxx[IJK(ix-2,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (sxy[IJK(ix,iy,iz)] - sxy[IJK(ix,iy-1,iz)]) - \
		  C0 * (sxy[IJK(ix,iy+1,iz)] - sxy[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();
	
	df_dK = ( C1 * (sxz[IJK(ix,iy,iz)] - sxz[IJK(ix,iy,iz-1)]) - \
		  C0 * (sxz[IJK(ix,iy,iz+1)] - sxz[IJK(ix,iy,iz-2)]) ) / SDMGeom->Dz();

	if (PROPAGATION == 0) {

	dsxx_dx[IJK(ix,iy,iz)] = pml_x->B_HALF(Gindx.x) * dsxx_dx[IJK(ix,iy,iz)] \
	  + pml_x->A_HALF(Gindx.x) * df_dI;

	dsxy_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dsxy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dsxz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dsxz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dsxx_dx[IJK(ix,iy,iz)];
	df_dJ += dsxy_dy[IJK(ix,iy,iz)];
	df_dK += dsxz_dz[IJK(ix,iy,iz)];

	}


	if (PROPAGATION == 1){
	  ux[IJK(ix,iy,iz)] = ux[IJK(ix,iy,iz)] + sgn * vx[IJK(ix,iy,iz)] * dt;
	}
	
	vx[IJK(ix,iy,iz)] = vx[IJK(ix,iy,iz)] + sgn * (dt / rho[IJK(ix,iy,iz)]) * \
	  (df_dI + df_dJ + df_dK);

	if (PROPAGATION == 0){
	  ux[IJK(ix,iy,iz)] = ux[IJK(ix,iy,iz)] + sgn * vx[IJK(ix,iy,iz)] * dt;

	}


	
      }
    }
  }


  

  
}


void SDM::Free_VY(VecI Init,VecI Iend, int zh){
  int iz = Iend.z + HALO.z - zh;
  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat rho_avg;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,rho_avg,Lindx,Gindx)\
  firstprivate(zh,iz)
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxy[IJK(ix+1,iy,iz)] - sxy[IJK(ix,iy,iz)]) - \
		  C0 * (sxy[IJK(ix+2,iy,iz)] - sxy[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (syy[IJK(ix,iy+1,iz)] - syy[IJK(ix,iy,iz)]) - \
		  C0 * (syy[IJK(ix,iy+2,iz)] - syy[IJK(ix,iy-1,iz)]) ) / SDMGeom->Dy();


	if (zh == 1) { // Free Surface Z = 0;
	
	  df_dK = -( (35.0 / 8.0) * syz[IJK(ix,iy,iz-1)] \
		    - (35.0 / 24.0) * syz[IJK(ix,iy,iz-2)] \
		    + (21.0 / 40.0) * syz[IJK(ix,iy,iz-3)] \
		    - (5.0 / 56.0) * syz[IJK(ix,iy,iz-4)] ) / SDMGeom->Dz();

	} else if (zh == 2) { // Free Surface Z == h;


	  df_dK = -( - (31.0 / 24.0) * syz[IJK(ix,iy,iz)] \
		    + (29.0 / 24.0) * syz[IJK(ix,iy,iz-1)] \
		    - (3.0 / 40.0) * syz[IJK(ix,iy,iz-2)] \
		    + (1.0 / 168.0) * syz[IJK(ix,iy,iz-3)] ) / SDMGeom->Dz();
	  
	}	


	if (PROPAGATION == 0){

	dsxy_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dsxy_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dsyy_dy[IJK(ix,iy,iz)] = pml_y->B(Gindx.y) * dsyy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A(Gindx.y) * df_dJ;

	dsyz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dsyz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dsxy_dx[IJK(ix,iy,iz)];
	df_dJ += dsyy_dy[IJK(ix,iy,iz)];
	df_dK += dsyz_dz[IJK(ix,iy,iz)];

	}
	
	rho_avg = ( rho[IJK(ix,iy,iz)] + rho[IJK(ix+1,iy,iz)] + \
		rho[IJK(ix,iy+1,iz)] + rho[IJK(ix+1,iy+1,iz)]) / 4.0;


	if (PROPAGATION == 1){
	  uy[IJK(ix,iy,iz)] = uy[IJK(ix,iy,iz)] + sgn * vy[IJK(ix,iy,iz)] * dt;
	}

	vy[IJK(ix,iy,iz)] = vy[IJK(ix,iy,iz)] + sgn *  (dt / rho_avg) * \
	  (df_dI + df_dJ + df_dK);


	if (PROPAGATION == 0){
	  uy[IJK(ix,iy,iz)] = uy[IJK(ix,iy,iz)] + sgn * vy[IJK(ix,iy,iz)] * dt;

	}



	
      }
    }

  
}


void SDM::FD_VY(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat rho_avg;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,rho_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxy[IJK(ix+1,iy,iz)] - sxy[IJK(ix,iy,iz)]) - \
		  C0 * (sxy[IJK(ix+2,iy,iz)] - sxy[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (syy[IJK(ix,iy+1,iz)] - syy[IJK(ix,iy,iz)]) - \
		  C0 * (syy[IJK(ix,iy+2,iz)] - syy[IJK(ix,iy-1,iz)]) ) / SDMGeom->Dy();
	
	df_dK = ( C1 * (syz[IJK(ix,iy,iz)] - syz[IJK(ix,iy,iz-1)]) - \
		  C0 * (syz[IJK(ix,iy,iz+1)] - syz[IJK(ix,iy,iz-2)]) ) / SDMGeom->Dz();


	if (PROPAGATION == 0) {

	dsxy_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dsxy_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dsyy_dy[IJK(ix,iy,iz)] = pml_y->B(Gindx.y) * dsyy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A(Gindx.y) * df_dJ;

	dsyz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dsyz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dsxy_dx[IJK(ix,iy,iz)];
	df_dJ += dsyy_dy[IJK(ix,iy,iz)];
	df_dK += dsyz_dz[IJK(ix,iy,iz)];

	}

	rho_avg = ( rho[IJK(ix,iy,iz)] + rho[IJK(ix+1,iy,iz)] + \
		rho[IJK(ix,iy+1,iz)] + rho[IJK(ix+1,iy+1,iz)]) / 4.0;



	if (PROPAGATION == 1){
	  uy[IJK(ix,iy,iz)] = uy[IJK(ix,iy,iz)] + sgn * vy[IJK(ix,iy,iz)] * dt;
	}

	vy[IJK(ix,iy,iz)] = vy[IJK(ix,iy,iz)] + sgn * (dt / rho_avg) * \
	  (df_dI + df_dJ + df_dK);


	if (PROPAGATION == 0){
	  uy[IJK(ix,iy,iz)] = uy[IJK(ix,iy,iz)] + sgn * vy[IJK(ix,iy,iz)] * dt;

	}



	
      }
    }
  }


  

  
}

void SDM::Free_VZ(VecI Init,VecI Iend){

  int iz = Iend.z + HALO.z;
  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat rho_avg;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,rho_avg,Lindx,Gindx)\
  firstprivate(iz)
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxz[IJK(ix+1,iy,iz)] - sxz[IJK(ix,iy,iz)]) - \
		  C0 * (sxz[IJK(ix+2,iy,iz)] - sxz[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (syz[IJK(ix,iy,iz)] - syz[IJK(ix,iy-1,iz)]) - \
		  C0 * (syz[IJK(ix,iy+1,iz)] - syz[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();


	df_dK = -( - (11.0 / 12.0) * szz[IJK(ix,iy,iz+1)] \
		  + (17.0 / 24.0) * szz[IJK(ix,iy,iz)] \
		  + (3.0 / 8.0) * szz[IJK(ix,iy,iz-1)]\
		  - (5.0 / 24.0) * szz[IJK(ix,iy,iz-2)] \
		  + (1.0 / 24.0) * szz[IJK(ix,iy,iz-3)] ) / SDMGeom->Dz();
	
	

	if (PROPAGATION == 0) {

	dsxz_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dsxz_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dsyz_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dsyz_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dszz_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dszz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A(Gindx.z) * df_dK;

	df_dI += dsxz_dx[IJK(ix,iy,iz)];
	df_dJ += dsyz_dy[IJK(ix,iy,iz)];
	df_dK += dszz_dz[IJK(ix,iy,iz)];

	}

	rho_avg = ( rho[IJK(ix,iy,iz)] + rho[IJK(ix+1,iy,iz)] + \
		rho[IJK(ix,iy,iz+1)] + rho[IJK(ix+1,iy,iz+1)]) / 4.0;


	if (PROPAGATION == 1){
	  uz[IJK(ix,iy,iz)] = uz[IJK(ix,iy,iz)] + sgn * vz[IJK(ix,iy,iz)] * dt;
	}


	vz[IJK(ix,iy,iz)] = vz[IJK(ix,iy,iz)] + sgn * (dt / rho_avg) * \
	  (df_dI + df_dJ + df_dK);


	if (PROPAGATION == 0){
	  uz[IJK(ix,iy,iz)] = uz[IJK(ix,iy,iz)] + sgn * vz[IJK(ix,iy,iz)] * dt;

	}


	
      }
    }
  
}


void SDM::FD_VZ(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat rho_avg;
  

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,rho_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (sxz[IJK(ix+1,iy,iz)] - sxz[IJK(ix,iy,iz)]) - \
		  C0 * (sxz[IJK(ix+2,iy,iz)] - sxz[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (syz[IJK(ix,iy,iz)] - syz[IJK(ix,iy-1,iz)]) - \
		  C0 * (syz[IJK(ix,iy+1,iz)] - syz[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();
	
	df_dK = ( C1 * (szz[IJK(ix,iy,iz+1)] - szz[IJK(ix,iy,iz)]) - \
		  C0 * (szz[IJK(ix,iy,iz+2)] - szz[IJK(ix,iy,iz-1)]) ) / SDMGeom->Dz();


	if (PROPAGATION == 0) {
	  
	dsxz_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dsxz_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dsyz_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dsyz_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dszz_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dszz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A(Gindx.z) * df_dK;

	df_dI += dsxz_dx[IJK(ix,iy,iz)];
	df_dJ += dsyz_dy[IJK(ix,iy,iz)];
	df_dK += dszz_dz[IJK(ix,iy,iz)];

	}
	
	rho_avg = ( rho[IJK(ix,iy,iz)] + rho[IJK(ix+1,iy,iz)] + \
		rho[IJK(ix,iy,iz+1)] + rho[IJK(ix+1,iy,iz+1)]) / 4.0;



	if (PROPAGATION == 1){
	  uz[IJK(ix,iy,iz)] = uz[IJK(ix,iy,iz)] + sgn * vz[IJK(ix,iy,iz)] * dt;
	}

	vz[IJK(ix,iy,iz)] = vz[IJK(ix,iy,iz)] + sgn *  (dt / rho_avg) * \
	  (df_dI + df_dJ + df_dK);


	if (PROPAGATION == 0){
	  uz[IJK(ix,iy,iz)] = uz[IJK(ix,iy,iz)] + sgn * vz[IJK(ix,iy,iz)] * dt;
	}


	
	
      }
    }
  }

  
}



void SDM::FD_SXY(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg;

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dJ = ( C1 * (vx[IJK(ix,iy+1,iz)] - vx[IJK(ix,iy,iz)]) - \
		  C0 * (vx[IJK(ix,iy+2,iz)] - vx[IJK(ix,iy-1,iz)]) ) / SDMGeom->Dy();
	
	df_dI = ( C1 * (vy[IJK(ix,iy,iz)] - vy[IJK(ix-1,iy,iz)]) - \
		  C0 * (vy[IJK(ix+1,iy,iz)] - vy[IJK(ix-2,iy,iz)]) ) / SDMGeom->Dx();
	

	if (PROPAGATION == 0) {
	  
	dvy_dx[IJK(ix,iy,iz)] = pml_x->B_HALF(Gindx.x) * dvy_dx[IJK(ix,iy,iz)] \
	  + pml_x->A_HALF(Gindx.x) * df_dI;

	dvx_dy[IJK(ix,iy,iz)] = pml_y->B(Gindx.y) * dvx_dy[IJK(ix,iy,iz)] \
	  + pml_y->A(Gindx.y) * df_dJ;

	df_dI += dvy_dx[IJK(ix,iy,iz)];
	df_dJ += dvx_dy[IJK(ix,iy,iz)];


	}
      
	mu_avg = 1.0 / ((0.50 / mu[IJK(ix,iy,iz)]) + (0.50 / mu[IJK(ix,iy+1,iz)]));


	sxy[IJK(ix,iy,iz)] = sxy[IJK(ix,iy,iz)] + sgn * (dt * mu_avg) * \
	  (df_dI + df_dJ);


	
      }
    }
  }

}


void SDM::Free_SXZ(VecI Init,VecI Iend){

  int iz = Iend.z + HALO.z;
  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg;

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,Lindx,Gindx)\
  firstprivate(iz)
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (vz[IJK(ix,iy,iz)] - vz[IJK(ix-1,iy,iz)]) - \
		  C0 * (vz[IJK(ix+1,iy,iz)] - vz[IJK(ix-2,iy,iz)]) ) / SDMGeom->Dx();


	df_dK = -( - (11.0 / 12.0) * vx[IJK(ix,iy,iz+1)] \
		  + (17.0 / 24.0) * vx[IJK(ix,iy,iz)] \
		  + (3.0 / 8.0) * vx[IJK(ix,iy,iz-1)]\
		  - (5.0 / 24.0) * vx[IJK(ix,iy,iz-2)] \
		  + (1.0 / 24.0) * vx[IJK(ix,iy,iz-3)] ) / SDMGeom->Dz();
	
	

	if (PROPAGATION == 0) {

	  dvz_dx[IJK(ix,iy,iz)] = pml_x->B_HALF(Gindx.x) * dvz_dx[IJK(ix,iy,iz)] \
	  + pml_x->A_HALF(Gindx.x) * df_dI;

	  dvx_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dvx_dz[IJK(ix,iy,iz)] \
	  + pml_z->A(Gindx.z) * df_dK;

	  df_dI += dvz_dx[IJK(ix,iy,iz)];
	  df_dK += dvx_dz[IJK(ix,iy,iz)];

	}

	
	mu_avg = 1.0 / ((0.50 / mu[IJK(ix,iy,iz)]) + (0.50 / mu[IJK(ix,iy,iz+1)]));

	sxz[IJK(ix,iy,iz)] = sxz[IJK(ix,iy,iz)] + sgn * (dt * mu_avg) * \
	  (df_dI + df_dK);


	
      }
    }

}

void SDM::FD_SXZ(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg;

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dI = ( C1 * (vz[IJK(ix,iy,iz)] - vz[IJK(ix-1,iy,iz)]) - \
		  C0 * (vz[IJK(ix+1,iy,iz)] - vz[IJK(ix-2,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dK = ( C1 * (vx[IJK(ix,iy,iz+1)] - vx[IJK(ix,iy,iz)]) - \
		  C0 * (vx[IJK(ix,iy,iz+2)] - vx[IJK(ix,iy,iz-1)]) ) / SDMGeom->Dz();
	


	if (PROPAGATION == 0) {
	  
	dvz_dx[IJK(ix,iy,iz)] = pml_x->B_HALF(Gindx.x) * dvz_dx[IJK(ix,iy,iz)] \
	  + pml_x->A_HALF(Gindx.x) * df_dI;

	dvx_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dvx_dz[IJK(ix,iy,iz)] \
	  + pml_z->A(Gindx.z) * df_dK;

	df_dI += dvz_dx[IJK(ix,iy,iz)];
	df_dK += dvx_dz[IJK(ix,iy,iz)];

	}

	mu_avg = 1.0 / ((0.50 / mu[IJK(ix,iy,iz)]) + (0.50 / mu[IJK(ix,iy,iz+1)]));

	sxz[IJK(ix,iy,iz)] = sxz[IJK(ix,iy,iz)] + sgn * (dt * mu_avg) * \
	  (df_dI + df_dK);


	
      }
    }
  }



  

  
}


void SDM::Free_SYZ(VecI Init,VecI Iend){

  int iz = Iend.z + HALO.z;
  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg;

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,Lindx,Gindx)\
  firstprivate(iz)
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dJ = ( C1 * (vz[IJK(ix,iy+1,iz)] - vz[IJK(ix,iy,iz)]) - \
		  C0 * (vz[IJK(ix,iy+2,iz)] - vz[IJK(ix,iy-1,iz)]) ) / SDMGeom->Dy();

	df_dK = -( - (11.0 / 12.0) * vy[IJK(ix,iy,iz+1)] \
		  + (17.0 / 24.0) * vy[IJK(ix,iy,iz)] \
		  + (3.0 / 8.0) * vy[IJK(ix,iy,iz-1)]\
		  - (5.0 / 24.0) * vy[IJK(ix,iy,iz-2)] \
		  + (1.0 / 24.0) * vy[IJK(ix,iy,iz-3)] ) / SDMGeom->Dz();
	


	if (PROPAGATION == 0) {
	  
	dvz_dy[IJK(ix,iy,iz)] = pml_y->B(Gindx.y) * dvz_dy[IJK(ix,iy,iz)] \
	  + pml_y->A(Gindx.y) * df_dJ;

	dvy_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dvy_dz[IJK(ix,iy,iz)] \
	  + pml_z->A(Gindx.z) * df_dK;

	df_dJ += dvz_dy[IJK(ix,iy,iz)];
	df_dK += dvy_dz[IJK(ix,iy,iz)];

	}

	
	mu_avg = 1.0 / ((0.1250 / mu[IJK(ix,iy,iz)]) + (0.1250 / mu[IJK(ix+1,iy,iz)]) + \
		(0.1250 / mu[IJK(ix,iy+1,iz)]) + (0.1250 / mu[IJK(ix+1,iy+1,iz)]) + \
		(0.1250 / mu[IJK(ix,iy,iz+1)]) + (0.1250 / mu[IJK(ix+1,iy,iz+1)]) + \
		(0.1250 / mu[IJK(ix,iy+1,iz+1)]) + (0.1250 / mu[IJK(ix+1,iy+1,iz+1)]));


	

	syz[IJK(ix,iy,iz)] = syz[IJK(ix,iy,iz)] + sgn * (dt * mu_avg) * \
	  (df_dJ + df_dK);


	
      }
    }

  
}


void SDM::FD_SYZ(VecI Init,VecI Iend){

  VecI Lindx,Gindx;
  Dfloat df_dI,df_dJ,df_dK;
  Dfloat mu_avg;

#pragma omp parallel for num_threads(N_omp)\
  private(df_dI,df_dJ,df_dK,mu_avg,Lindx,Gindx)
  for (int iz=Init.z + HALO.z; iz<Iend.z + HALO.z; ++iz){
    for (int iy=Init.y + HALO.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=Init.x + HALO.x; ix<Iend.x + HALO.x; ++ix){


	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);

	df_dJ = ( C1 * (vz[IJK(ix,iy+1,iz)] - vz[IJK(ix,iy,iz)]) - \
		  C0 * (vz[IJK(ix,iy+2,iz)] - vz[IJK(ix,iy-1,iz)]) ) / SDMGeom->Dy();
	
	df_dK = ( C1 * (vy[IJK(ix,iy,iz+1)] - vy[IJK(ix,iy,iz)]) - \
		  C0 * (vy[IJK(ix,iy,iz+2)] - vy[IJK(ix,iy,iz-1)]) ) / SDMGeom->Dz();
	


	if (PROPAGATION == 0) {
	  
	dvz_dy[IJK(ix,iy,iz)] = pml_y->B(Gindx.y) * dvz_dy[IJK(ix,iy,iz)] \
	  + pml_y->A(Gindx.y) * df_dJ;

	dvy_dz[IJK(ix,iy,iz)] = pml_z->B(Gindx.z) * dvy_dz[IJK(ix,iy,iz)] \
	  + pml_z->A(Gindx.z) * df_dK;

	df_dJ += dvz_dy[IJK(ix,iy,iz)];
	df_dK += dvy_dz[IJK(ix,iy,iz)];

	}
	
	
	mu_avg = 1.0 / ((0.1250 / mu[IJK(ix,iy,iz)]) + (0.1250 / mu[IJK(ix+1,iy,iz)]) + \
                (0.1250 / mu[IJK(ix,iy+1,iz)]) + (0.1250 / mu[IJK(ix+1,iy+1,iz)]) + \
                (0.1250 / mu[IJK(ix,iy,iz+1)]) + (0.1250 / mu[IJK(ix+1,iy,iz+1)]) + \
                (0.1250 / mu[IJK(ix,iy+1,iz+1)]) + (0.1250 / mu[IJK(ix+1,iy+1,iz+1)]));


	syz[IJK(ix,iy,iz)] = syz[IJK(ix,iy,iz)] + sgn *  (dt * mu_avg) * \
	  (df_dJ + df_dK);


	
      }
    }
  }

  
}


void SDM::Free_SII(VecI Init,VecI Iend, int zh){

  int iz = Iend.z + HALO.z - zh ;
  Dfloat d_free;
  VecI Lindx,Gindx;
  Dfloat mu_avg,lamb_avg,df_dI,df_dJ,df_dK;
  Dfloat df_dI_free,df_dJ_free;

   // Free Surface Implementation Stress Imaging
  
  #pragma omp parallel for num_threads(N_omp)\
    private(df_dI,df_dJ,df_dI_free,df_dJ_free,df_dK,mu_avg,lamb_avg,Lindx,Gindx,d_free)	\
    firstprivate(iz,zh)
   for (int iy=HALO.y+Init.y; iy<Iend.y + HALO.y; ++iy){
      for (int ix=HALO.x+Init.x; ix<Iend.x + HALO.x; ++ix){

	Lindx = {ix,iy,iz};
	Gindx = Loc2Glo(Lindx);


	// Average properties

	mu_avg = 1.0 / ((0.50 / mu[IJK(ix,iy,iz)]) + (0.50 / mu[IJK(ix+1,iy,iz)]));
	lamb_avg = 1.0 / ((0.50 / lamb[IJK(ix,iy,iz)]) + (0.50 / lamb[IJK(ix+1,iy,iz)]));

	

	df_dI = ( C1 * (vx[IJK(ix+1,iy,iz)] - vx[IJK(ix,iy,iz)]) - \
		  C0 * (vx[IJK(ix+2,iy,iz)] - vx[IJK(ix-1,iy,iz)]) ) / SDMGeom->Dx();
	
	df_dJ = ( C1 * (vy[IJK(ix,iy,iz)] - vy[IJK(ix,iy-1,iz)]) - \
		  C0 * (vy[IJK(ix,iy+1,iz)] - vy[IJK(ix,iy-2,iz)]) ) / SDMGeom->Dy();

	
	if (zh == 1) { // z = 0; Free Surface
	
	df_dK = - (df_dI + df_dJ) * (lamb_avg / (lamb_avg + 2.0 * mu_avg));

	} else if (zh == 2) { // z = h; Free Surface


	df_dI_free = ( C1 * (vx[IJK(ix+1,iy,iz+1)] - vx[IJK(ix,iy,iz+1)]) - \
		  C0 * (vx[IJK(ix+2,iy,iz+1)] - vx[IJK(ix-1,iy,iz+1)]) ) / SDMGeom->Dx();
	
	df_dJ_free = ( C1 * (vy[IJK(ix,iy,iz+1)] - vy[IJK(ix,iy-1,iz+1)]) - \
		  C0 * (vy[IJK(ix,iy+1,iz+1)] - vy[IJK(ix,iy-2,iz+1)]) ) / SDMGeom->Dy();


	  d_free = - (df_dI_free + df_dJ_free) * (lamb_avg / (lamb_avg + 2.0 * mu_avg)) * SDMGeom->Dz();

	  df_dK = -(1.0 / SDMGeom->Dz()) * ( - (1.0 / 22.0) * d_free - (577.0 / 528.0) * vz[IJK(ix,iy,iz)] \
				    + (201.0 / 176.0) * vz[IJK(ix,iy,iz-1)] \
				    - (9.0 / 176.0) * vz[IJK(ix,iy,iz-2)] \
					    + (1.0 / 528.0) * vz[IJK(ix,iy,iz-3)] );
	}



	if (PROPAGATION == 0) {
	  
	dvx_dx[IJK(ix,iy,iz)] = pml_x->B(Gindx.x) * dvx_dx[IJK(ix,iy,iz)] \
	  + pml_x->A(Gindx.x) * df_dI;

	dvy_dy[IJK(ix,iy,iz)] = pml_y->B_HALF(Gindx.y) * dvy_dy[IJK(ix,iy,iz)] \
	  + pml_y->A_HALF(Gindx.y) * df_dJ;

	dvz_dz[IJK(ix,iy,iz)] = pml_z->B_HALF(Gindx.z) * dvz_dz[IJK(ix,iy,iz)]  \
	  + pml_z->A_HALF(Gindx.z) * df_dK;

	df_dI += dvx_dx[IJK(ix,iy,iz)];
	df_dJ += dvy_dy[IJK(ix,iy,iz)];
	df_dK += dvz_dz[IJK(ix,iy,iz)];

	}
	
	// SXX

	sxx[IJK(ix,iy,iz)] = sxx[IJK(ix,iy,iz)] +  sgn * \
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dI) + sgn * \
	  dt * lamb_avg *  (df_dJ + df_dK);

	// SYY
	
	syy[IJK(ix,iy,iz)] = syy[IJK(ix,iy,iz)] + sgn * \
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dJ) + sgn * \
	  dt * lamb_avg *  (df_dI + df_dK);


	// SZZ

	if (zh == 1) {
	  
	szz[IJK(ix,iy,iz)] = 0.0;
	
	} else if (zh == 2) {

	  szz[IJK(ix,iy,iz)] = szz[IJK(ix,iy,iz)] + sgn * \
	  (dt * (lamb_avg + 2.0 * mu_avg) * df_dK) + sgn * \
	  dt * lamb_avg *  (df_dJ + df_dI);
	}

	

      }
   }



}





void SDM::FDSII() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};

  
  if ((PROPAGATION == 0)||(PROPAGATION == -1)) {

  if (Nsdm.x == NumSubDom.x-1) { 
    end.x = NodLoc.x - 1;
  }

 
  if (Nsdm.z == NumSubDom.z-1){

    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_SII(init,end,1); // z = 0;
    Free_SII(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
    
  } 


  FD_SII(init,end);

  } else if (PROPAGATION == 1) {


    if (Nsdm.x == 0) {
      init.x = PML.x-1;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x + 1;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y + 1;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z;
    }

    if (Nsdm.z == NumSubDom.z-1){

    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_SII(init,end,1); // z = 0;
    Free_SII(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
    

    } 
    
    FD_SII(init,end);
    
  }
  
  
}


void SDM::FDSXY() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};

  if ((PROPAGATION == 0)||(PROPAGATION == -1)){
    
  if (Nsdm.y == NumSubDom.y - 1) { 
    end.y = NodLoc.y - 1;
  }
  
  FD_SXY(init,end);

  } else if (PROPAGATION == 1) {


   if (Nsdm.x == 0) {
      init.x = PML.x;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x + 1;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y - 1;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z;
    }
    
    FD_SXY(init,end);

  }

}

void SDM::FDSXZ() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};


  if ((PROPAGATION == 0)||(PROPAGATION == -1)){
    
  if (Nsdm.z == NumSubDom.z - 1) {
    // FREE SURFACE
    // H - AFDA Kristek(2002);
    end.z = NodLoc.z - 2;
    Free_SXZ(init,end);
  }
    
  FD_SXZ(init,end);

  } else if (PROPAGATION == 1){


    if (Nsdm.x == 0) {
      init.x = PML.x;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x + 1;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z - 1;
    }
    

    if (Nsdm.z == NumSubDom.z - 1) {
    // FREE SURFACE
    // H - AFDA Kristek(2002);
    end.z = NodLoc.z - 2;
    Free_SXZ(init,end);
  }
    
  FD_SXZ(init,end);

  }

}

void SDM::FDSYZ() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};


  if ((PROPAGATION == 0)||(PROPAGATION == -1)) {
    
   if (Nsdm.x == NumSubDom.x - 1) { 
    end.x = NodLoc.x - 1;
  }

   if (Nsdm.y == NumSubDom.y - 1) { 
    end.y = NodLoc.y - 1;
  }


   if (Nsdm.z == NumSubDom.z - 1) {

     // FREE SURFACE
    // H - AFDA Kristek(2002);
    end.z = NodLoc.z - 2;
    Free_SYZ(init,end);
    
   }
    
  FD_SYZ(init,end);

  } else if (PROPAGATION == 1) {

    if (Nsdm.x == 0) {
      init.x = PML.x;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y - 1;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z - 1;
    }


    if (Nsdm.z == NumSubDom.z - 1) {

     // FREE SURFACE
    // H - AFDA Kristek(2002);
    end.z = NodLoc.z - 2;
    Free_SYZ(init,end);
    
   }
    
  FD_SYZ(init,end);
    
    
  }

}

void SDM::FDVX() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};


  if ((PROPAGATION == 0)||(PROPAGATION == -1)){
    
  if (Nsdm.z == 0) { 
    init.z = 1;
  }

  if (Nsdm.x == 0) { 
    init.x = 1;
  }

  if (Nsdm.x == NumSubDom.x - 1) { 
    end.x = NodLoc.x - 1;
  }

  if (Nsdm.y == 0) { 
    init.y = 1;
  }

  if (Nsdm.y == NumSubDom.y - 1) { 
    end.y = NodLoc.y - 1;
  }

  
  if (Nsdm.z == NumSubDom.z-1){
    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_VX(init,end,1); // z = 0;
    Free_VX(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
  }

  FD_VX(init,end);

  } else if (PROPAGATION == 1) {

    if (Nsdm.x == 0) {
      init.x = PML.x;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x + 1;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y + 1;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z;
    }
    
    if (Nsdm.z == NumSubDom.z-1){
    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_VX(init,end,1); // z = 0;
    Free_VX(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
  }

  FD_VX(init,end);


  }

}

void SDM::FDVY() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};

  if ((PROPAGATION == 0)||(PROPAGATION == -1)){
    
  if (Nsdm.z == 0) { 
    init.z = 1;
  }

  if (Nsdm.y == NumSubDom.y - 1) { 
    end.y = NodLoc.y - 1;
  }

  if (Nsdm.x == NumSubDom.x - 1) { 
    end.x = NodLoc.x - 1;
  }


  if (Nsdm.z == NumSubDom.z-1){
    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_VY(init,end,1); // z = 0;
    Free_VY(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
  }
  
  FD_VY(init,end);

  } else if (PROPAGATION == 1) {


    if (Nsdm.x == 0) {
      init.x = PML.x - 1;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y - 1;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z;
    }


    if (Nsdm.z == NumSubDom.z-1){
    // FREE SURFACE 

    // H - AFDA Kristek(2002);
    Free_VY(init,end,1); // z = 0;
    Free_VY(init,end,2); // z = h;
    end.z = NodLoc.z - 2;
  }
  
  FD_VY(init,end);

  }

}

void SDM::FDVZ() {
  VecI init,end;

  init = {0,0,0};
  end = {NodLoc.x,NodLoc.y,NodLoc.z};


  if ((PROPAGATION == 0)||(PROPAGATION == -1)){
    
  if (Nsdm.y == NumSubDom.y - 1) { 
    end.y = NodLoc.y - 1;
  }

   if (Nsdm.y == 0) { 
    init.y = 1;
  }

    if (Nsdm.x == NumSubDom.x - 1) { 
    end.x = NodLoc.x - 1;
  }


    // FREE SURFACE H - AFDA KRISTEK(2002);
    if (Nsdm.z == NumSubDom.z - 1) { 
    end.z = NodLoc.z - 2;
    Free_VZ(init,end);

    //for (int iy= HALO.y; iy<SNodeY(); ++iy){
    //for (int ix=HALO.x; ix<SNodeX(); ++ix){
    //vz[IJK(ix,iy,SNodeZ()-KHALO-1)] = vz[IJK(ix,iy,SNodeZ()-KHALO-2)] ;
    //uz[IJK(ix,iy,SNodeZ()-KHALO-1)] = uz[IJK(ix,iy,SNodeZ()-KHALO-2)] ;
    // }
    //}
    
    }

    FD_VZ(init,end);

  } else if (PROPAGATION == 1) {


    if (Nsdm.x == 0) {
      init.x = PML.x - 1;
    }

    if (Nsdm.x == NumSubDom.x-1) {
      end.x = SNodeX() - PML.x - 2 * HALO.x;
    }

    if (Nsdm.y == 0) {
      init.y = PML.y;
    }

    if (Nsdm.y == NumSubDom.y-1) {
      end.y = SNodeY() - PML.y - 2 * HALO.y + 1;
    }

    if (Nsdm.z == 0) {
      init.z = PML.z - 1;
    }


    // FREE SURFACE H - AFDA KRISTEK(2002);
    if (Nsdm.z == NumSubDom.z - 1) { 
    end.z = NodLoc.z - 2;
    Free_VZ(init,end);
    }

    //for (int iy= HALO.y; iy<SNodeY(); ++iy){
    //for (int ix=HALO.x; ix<SNodeX(); ++ix){
    //vz[IJK(ix,iy,SNodeZ()-KHALO-1)] = vz[IJK(ix,iy,SNodeZ()-KHALO-2)] ;
    //uz[IJK(ix,iy,SNodeZ()-KHALO-1)] = uz[IJK(ix,iy,SNodeZ()-KHALO-2)] ;
    // }
    // }
    
    FD_VZ(init,end);

  }

}


void SDM::InitSource(geometry3D *GDomain,std::string nFile,int nsource,int SrcFile,int nt){

  FileSrcB = SrcFile;
  ntimesrc = nt; // Variable Write File Source
  srct = new Dfloat[nt];
  
  if (FileSrcB == 0) {
  sourceM = new source(GDomain,nFile,nsource);
  } else if (FileSrcB == 1){
    FileSrc = "SrcTime.bin";
    std::fstream R;
    R.open("SrcTime.bin",std::ios::binary | std::ios::in);
    for (int it = 0; it<nt;++it){
      R.read( (char*)&srct[it], sizeof(Dfloat));
    }

    R.close();
    sourceM = new source(GDomain,nFile,nsource);

  } else {
    printf("This is not a option for the source/n");
  }


  idx_source = new VecI[nsource];
	  
  for (int i = 0; i<sourceM->ns; ++i){
    idx_source[i] = SFindNode(sourceM->pos_src[i]);
  }
  
}


void SDM::InitRecept(geometry3D *GDomain,std::string nFile,int nrecep,int nt){
  
  station  = new receptor(GDomain,nFile,nrecep,nt);

  Rvx = new Dfloat[station->nt * station->nr];
  Rvy = new Dfloat[station->nt * station->nr];
  Rvz = new Dfloat[station->nt * station->nr];

  idx_station = new VecI[nrecep];

  for (int i = 0; i<station->nr; ++i){
    idx_station[i] = SFindNode(station->pos_recep[i]);
  }
  
 
}


void SDM::InitAdj(geometry3D *GDomain,std::string nFile,int nrecep,int nt){
  
  station  = new receptor(GDomain,nFile,nrecep,nt);

  idx_station = new VecI[nrecep];

  for (int i = 0; i<station->nr; ++i){
    idx_station[i] = SFindNode(station->pos_recep[i]);

      if ( (idx_station[i].x > -1) && (idx_station[i].y > -1) && \
	   (idx_station[i].z > -1) ){

	station->FileOpen(i,1);
	station->LoadFile(i);

      }
    
  }
  
 
}

void SDM::AddSourceAdj(int itime){
  Dfloat Rvx,Rvy,Rvz;
  Dfloat rho_cte;

   for (int i = 0; i<station->nr; ++i){
    
     
      if ( (idx_station[i].x > -1) && (idx_station[i].y > -1)  &&\
	   (idx_station[i].z > -1) ){


	rho_cte = 2700.0;

	AddVal(station->pos_recep[i],"VX", (dt / rho_cte) * station->vx_ad[itime + station->nt * i]);
	AddVal(station->pos_recep[i],"VY", (dt / rho_cte) * station->vy_ad[itime + station->nt * i]);
	AddVal(station->pos_recep[i],"VZ", (dt / rho_cte) * station->vz_ad[itime + station->nt * i]);
        
      }
   }
  

}




void SDM::EndRecept(){
  
  for (int i = 0; i<station->nr; ++i){
 
      if ( (idx_station[i].x > -1) && (idx_station[i].y > -1) &&\
	   (idx_station[i].z > -1) ){

	station->FileOpen(i,PROPAGATION);

	//for (int ktime=0; ktime<station->nt; ++ktime){

	station->WriteFile(i,Rvx + i * station->nt,\
			   Rvy + i * station->nt,\
			   Rvz + i * station->nt);
	//}

	
	station->FileClose(i);
	

      }
    
  }

  delete [] Rvx;
  delete [] Rvy;
  delete [] Rvz;
  delete [] idx_station;
  delete station;
  
}


void SDM::GetRecept(int ktime){

   for (int i = 0; i<station->nr; ++i){
    
     
      if ( (idx_station[i].x > -1) && (idx_station[i].y > -1)  &&\
	   (idx_station[i].z > -1) ){

	Rvx[ktime + i * station->nt] = GetVal(station->pos_recep[i],"VX");
	Rvy[ktime + i * station->nt] = GetVal(station->pos_recep[i],"VY");
	Rvz[ktime + i * station->nt] = GetVal(station->pos_recep[i],"VZ");

	//station->WriteFile(i,Rvx,Rvy,Rvz);
	
      }
    
   }
  

}



void SDM::AddSource(int itime, int T_SRC){

  Dfloat st = 0;
  Dfloat st_sinc;


  for (int i = 0; i<sourceM->ns; ++i){

    if (FileSrcB == 0){
    st = sourceM->sourceType(sourceM->d_t0[i],f0,itime,dt,T_SRC);
    } else if (FileSrcB == 1){
      st = srct[itime];
    }

    // SAVE SOURCE TIME FUNCTION FILE 
    srct[itime] = st;
    
    VecI b = sourceM->pos_src[i];


    // Vertical Source
    if (sourceM->M0[i] < 0) {

      Dfloat rho_cte = 2700.0;
      st *= (dt / rho_cte);
      AddVal({b.x,b.y,b.z},"VX", st * sourceM->Myy[i]*sgn);
      AddVal({b.x,b.y,b.z},"VY", st * sourceM->Mxx[i]*sgn);
      AddVal({b.x,b.y,b.z},"VZ", st * sourceM->Mzz[i]*sgn);
	
      // Moment Tensor Source
    } else if (sourceM->M0[i] > 0) {

      st_sinc = st * -dt;
   
      st *= (-dt / (SDMGeom->Dx() * SDMGeom->Dy() * SDMGeom->Dz()) );

    // CLOSE FREE SURFACE SOURCE
    if (sourceM->nshift[i] < 4) {
      sourceM->w_sinc(1,i);
     // BELOW FREE SURFACE
    } else if (sourceM->nshift[i] >= 4){
      sourceM->w_sinc(0,i);
    }
   
// MOMENT TENSOR ORIENTATION
// X: EAST
// Y: NORTH
// Z: UPWARD

    // AddVal({b.x,b.y,b.z},"SXX", st * sourceM->Myy[i]*sgn);
    // AddVal({b.x,b.y,b.z},"SYY", st * sourceM->Mxx[i]*sgn);
    // AddVal({b.x,b.y,b.z},"SZZ", st * sourceM->Mzz[i]*sgn);
    // AddVal({b.x,b.y,b.z},"SXY", st * sourceM->Mxy[i]*sgn);
    // AddVal({b.x,b.y,b.z},"SXZ", -st * sourceM->Myz[i]*sgn);
    // AddVal({b.x,b.y,b.z},"SYZ", -st * sourceM->Mxz[i]*sgn);


    

    Dfloat w_sii,w_sxy,w_sxz;
    AddVal({b.x,b.y,b.z},"SYZ", -st * sourceM->Mxz[i]*sgn);
    
    for (int j = -3;j<=4;j++){
      for (int k = -3;k<=4;k++){
	
    	w_sii = sourceM->sinc_wy[3 + j] * sourceM->sinc_wz[3 + k] * (1.0 / SDMGeom->Dx());
    	w_sxy = sourceM->sinc_wx[3 + j] * sourceM->sinc_wz[3 + k] * (1.0 / SDMGeom->Dy());
    	w_sxz = sourceM->sinc_wx[3 + j] * sourceM->sinc_wy[3 + k] * (1.0 / SDMGeom->Dz());
	
    	AddVal({b.x,b.y+j,b.z+k},"SXX", st_sinc * sourceM->Myy[i]*sgn * w_sii);
    	AddVal({b.x,b.y+j,b.z+k},"SYY", st_sinc * sourceM->Mxx[i]*sgn * w_sii);
    	AddVal({b.x,b.y+j,b.z+k},"SZZ", st_sinc * sourceM->Mzz[i]*sgn * w_sii);
    	AddVal({b.x+j,b.y,b.z+k},"SXY", st_sinc * sourceM->Mxy[i]*sgn * w_sxy);
    	AddVal({b.x+j,b.y+k,b.z},"SXZ", -st_sinc * sourceM->Myz[i] * sgn * w_sxz);
      }
    }
    
  
    }

  }
  

}


void SDM::EndSource(){
  std::fstream R;

  R.open("SOURCE.bin",std::ios::binary | std::ios::out);

  for (int it = 0; it<ntimesrc;++it){
    R.write( (char*)&srct[it], sizeof(Dfloat));
    }

  R.close();

  delete [] srct;
  
}



void SDM::printfile(Dfloat *Var,char *nfile,int ktime){
  FILE *R;
  char times[200];
  MPI_File fhw;
  MPI_Status status;
  int offset;

  int subindx = Nsdm.x + Nsdm.y * NumSubDom.x + Nsdm.z * NumSubDom.x * NumSubDom.y;

  
  sprintf(times,"temp/%s-%d.bin",nfile,ktime);

  MPI_File_open(MPI_COMM_WORLD,times,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fhw);

  offset = subindx * SDMGeom->HALO_Node() * sizeof(Dfloat);

  MPI_File_write_at(fhw,offset,Var,SDMGeom->HALO_Node(),MY_MPI_Dfloat,&status);
  
  MPI_File_close(&fhw);

  
  
   // R=fopen(times,"wb");
   // fwrite(Var,sizeof(Dfloat),SDMGeom->HALO_Node(),R);
   
   // for (int iz=0;iz<SDMGeom->L_NodeZ();iz++){
   //  for (int iy=0;iy<SDMGeom->L_NodeY();iy++){
   //    for (int ix=0;ix<SDMGeom->L_NodeX();ix++){
   // int indx = (ix + HALO.x) + (iy + HALO.y) * SDMGeom->HALO_NodeX() +	\
   //  	   (iz + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();

   // 	 fwrite(&Var[indx],sizeof(Dfloat),1,R);


   //  }
   // }
   // }

    // fclose(R);	
    

  
}



void SDM::loadfile(Dfloat *Var,char *nfile,int ktime){
  // FILE *R;
  // char times[200];

  // int subindx = Nsdm.x + Nsdm.y * NumSubDom.x + Nsdm.z * NumSubDom.x * NumSubDom.y;
  // sprintf(times,"temp/%s_%d-%d.bin",nfile,subindx,ktime);
  
  
  //  R=fopen(times,"rb");

  //  fread(Var,sizeof(Dfloat),SDMGeom->HALO_Node(),R);
   
  //  // for (int iz=0;iz<SDMGeom->L_NodeZ();iz++){
  //  // for (int iy=0;iy<SDMGeom->L_NodeY();iy++){
  //  //     for (int ix=0;ix<SDMGeom->L_NodeX();ix++){
  //  // 	 int indx = (ix + HALO.x) + (iy + HALO.y) * SDMGeom->HALO_NodeX() + \
  //  //  	   (iz + HALO.z) * SDMGeom->HALO_NodeX() * SDMGeom->HALO_NodeY();

  //  // 	 fread(&Var[indx],sizeof(Dfloat),1,R);


  //  // 	 }
  //  //     }
  //  //   }

  //   fclose(R);



  FILE *R;
  char times[200];
  MPI_File fhw;
  MPI_Status status;
  int offset;

  int subindx = Nsdm.x + Nsdm.y * NumSubDom.x + Nsdm.z * NumSubDom.x * NumSubDom.y;

  
  sprintf(times,"temp/%s-%d.bin",nfile,ktime);

  MPI_File_open(MPI_COMM_WORLD,times,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&fhw);

  offset = subindx * SDMGeom->HALO_Node() * sizeof(Dfloat);

  MPI_File_read_at(fhw,offset,Var,SDMGeom->HALO_Node(),MY_MPI_Dfloat,&status);
  
  MPI_File_close(&fhw);

  
}

void SDM::file(char *NameVar, int time ,int oi){


   if (strcmp("SXX",NameVar) == 0){

     if (oi == 0)
       printfile(sxx,NameVar,time);

     if (oi == 1)
       loadfile(sxx,NameVar,time);
       
  }
  if (strcmp("SYY",NameVar) == 0){

    if (oi == 0)
      printfile(syy,NameVar,time);

    if (oi == 1)
      loadfile(syy,NameVar,time);
  }

  if (strcmp("SZZ",NameVar) == 0){

    if (oi == 0)
      printfile(szz,NameVar,time);

    if (oi == 1)
      loadfile(szz,NameVar,time);
  }

  if (strcmp("SXY",NameVar) == 0){

    if (oi == 0)
      printfile(sxy,NameVar,time);

    if (oi == 1)
      loadfile(sxy,NameVar,time);
  }

  if (strcmp("SXZ",NameVar) == 0){

    if (oi == 0)
    printfile(sxz,NameVar,time);

    if (oi == 1)
    loadfile(sxz,NameVar,time);
  }

  if (strcmp("SYZ",NameVar) == 0){

    if (oi == 0)
      printfile(syz,NameVar,time);

    if (oi == 1)
      loadfile(syz,NameVar,time);
  }

  if (strcmp("VX",NameVar) == 0){

    if (oi == 0)   
      printfile(vx,NameVar,time);

    if (oi == 1)   
      loadfile(vx,NameVar,time);
  }

  if (strcmp("VY",NameVar) == 0){

    if (oi == 0)
      printfile(vy,NameVar,time);

    if (oi == 1)
      loadfile(vy,NameVar,time);
  }

  if (strcmp("VZ",NameVar) == 0){

    if (oi == 0)
      printfile(vz,NameVar,time);

    if (oi == 1)
      loadfile(vz,NameVar,time);
  }


    if (strcmp("UX",NameVar) == 0){

    if (oi == 0)   
      printfile(ux,NameVar,time);

    if (oi == 1)   
      loadfile(ux,NameVar,time);
  }

  if (strcmp("UY",NameVar) == 0){

    if (oi == 0)
      printfile(uy,NameVar,time);

    if (oi == 1)
      loadfile(uy,NameVar,time);
  }

  if (strcmp("UZ",NameVar) == 0){

    if (oi == 0)
      printfile(uz,NameVar,time);

    if (oi == 1)
      loadfile(uz,NameVar,time);
  }


}

int SDM::BNorth(){

if (Nsdm.y == NumSubDom.y-1){
return 0;
}else {
  return 1;
} 
}

int SDM::BSouth(){

if (Nsdm.y == 0){
return 0;
}else {
  return 1;
} 
}

int SDM::BWest(){

if (Nsdm.x == 0){
return 0;
}else {
  return 1;
} 
}


int SDM::BEast(){

if (Nsdm.x == NumSubDom.x-1){
return 0;
}else {
  return 1;
} 
}

int SDM::BDown(){

if (Nsdm.z == 0){
return 0;
}else {
  return 1;
} 
}

int SDM::BUp(){

if (Nsdm.z == NumSubDom.z-1){
return 0;
}else {
  return 1;
} 
}

void SDM::boundX(Dfloat *var,int side,int inout,int time,int indx_var){

  int indx,indx1,size;
  size = KHALO * SDMGeom->L_NodeY() * SDMGeom->L_NodeZ();
  // Left Side

   // SAVE BOUNDARY
   if (inout == 0) {

  if (side == 0){
    
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<KHALO;i++){
	  indx = IJK(PML.x + i,HALO.y + j,HALO.z + k);
	  indx1 = i + j * KHALO + k * KHALO * SDMGeom->L_NodeY();
	  bn_lx[(indx_var - 1) * size + indx1 + time * size * 9] = var[indx];
	}
      }
    }

  

// Right Side
} else if (side == 1) {


    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<KHALO;i++){

	  
	  indx = IJK(SNodeX() - KHALO - PML.x + i,HALO.y + j,HALO.z + k);
	  indx1 = i + j * KHALO + k * KHALO * SDMGeom->L_NodeY();
	  bn_rx[(indx_var - 1) * size + indx1 + time * size * 9] = var[indx];
	  
	}
      }
    }
    
      }

  // READ BOUNDARY
   } else if (inout == 1) {



      if (side == 0){

     for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<KHALO;i++){
	  
	  indx = IJK(PML.x + i,HALO.y + j,HALO.z + k);
	  indx1 = i + j * KHALO + k * KHALO * SDMGeom->L_NodeY();
	  var[indx] = bn_lx[(indx_var - 1) * size + indx1 + time * size * 9];
	}
      }
    }
  

// Right Side
} else if (side == 1) {


    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<KHALO;i++){

	  indx = IJK(SNodeX() - KHALO - PML.x + i,HALO.y + j,HALO.z + k);
	  indx1 = i + j * KHALO + k * KHALO * SDMGeom->L_NodeY();
	  var[indx] = bn_rx[(indx_var - 1) * size + indx1 + time * size * 9];
	  
	}
      }
    }


    
      }

   }

}

void SDM::boundY(Dfloat *var,int side,int inout,int time,int indx_var){
  int indx,indx1,size;
  size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeZ();

   // SAVE BOUNDARY
   if (inout == 0) {

  if (side == 0){

    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<KHALO;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){
	  indx = IJK(HALO.x + i,PML.y + j,HALO.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k * KHALO * SDMGeom->L_NodeX();
	  bn_ly[(indx_var - 1) * size + indx1 + time * size * 9] = var[indx];

	  
	}
      }
    }
    

// Right Side
} else if (side == 1) {

  


    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<KHALO;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){
	  
	  indx = IJK(HALO.x + i,SNodeY() - KHALO - PML.y + j,HALO.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k * KHALO * SDMGeom->L_NodeX();
	  bn_ry[(indx_var - 1) * size + indx1 + time * size * 9] = var[indx];
	  
	}
      }
    }

  
    
      }

  // READ BOUNDARY
   } else if (inout == 1) {



      if (side == 0){
	
    
    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<KHALO;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){
	  
	  indx = IJK(HALO.x + i,PML.y + j,HALO.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k * KHALO * SDMGeom->L_NodeX();
	  var[indx] = bn_ly[(indx_var - 1) * size + indx1 + time * size * 9];
	  
	}
      }
    }


// Right Side
} else if (side == 1) {




    for (int k=0;k<SDMGeom->L_NodeZ();k++){
      for (int j=0;j<KHALO;j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){
	  
	  indx = IJK(HALO.x + i,SNodeY() - KHALO - PML.y + j,HALO.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k * KHALO * SDMGeom->L_NodeX();
	  var[indx] = bn_ry[(indx_var - 1) * size + indx1 + time * size * 9];
	  
	}
      }
    }
    
      }

   }

}


void SDM::boundZ(Dfloat *var,int side,int inout,int time,int indx_var){

  int indx,indx1,size;
  size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeY();
  // Left Side

   FILE *R;
   char times[200];


   // SAVE BOUNDARY
   if (inout == 0) {

  if (side == 0){

    
    for (int k=0;k<KHALO;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){
	  
	  indx = IJK(HALO.x + i,HALO.y + j,PML.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k *  SDMGeom->L_NodeY() * SDMGeom->L_NodeX();
	  bn_lz[(indx_var - 1) * size + indx1 + time * size * 9] = var[indx];
	  
	}
      }
    }


  }

  // READ BOUNDARY
   } else if (inout == 1) {



      if (side == 0){




    for (int k=0;k<KHALO;k++){
      for (int j=0;j<SDMGeom->L_NodeY();j++){
	for (int i=0;i<SDMGeom->L_NodeX();i++){

	  indx = IJK(HALO.x + i,HALO.y + j,PML.z + k);
	  indx1 = i + j * SDMGeom->L_NodeX() + k *  SDMGeom->L_NodeY() * SDMGeom->L_NodeX();
	  var[indx] = bn_lz[(indx_var - 1) * size + indx1 + time * size * 9];
	}
      }
    }

      } 

   }


}


void SDM::SaveBoundaries_V(int itime){

  int time;
  
      time = itime - nsteps * int(itime/nsteps);

      //printf("%d\n",time);


  if (Nsdm.x == 0){

    boundX(vx,0,0,time,1);
    boundX(vy,0,0,time,2);
    boundX(vz,0,0,time,3);

  }else if (Nsdm.x == NumSubDom.x-1) {

    boundX(vx,1,0,time,1);
    boundX(vy,1,0,time,2);
    boundX(vz,1,0,time,3);
  }

    if (Nsdm.y == 0){

    boundY(vx,0,0,time,1);
    boundY(vy,0,0,time,2);
    boundY(vz,0,0,time,3);

  }else if (Nsdm.y == NumSubDom.y-1) {

    boundY(vx,1,0,time,1);
    boundY(vy,1,0,time,2);
    boundY(vz,1,0,time,3);
  }

    if (Nsdm.z == 0){

    boundZ(vx,0,0,time,1);
    boundZ(vy,0,0,time,2);
    boundZ(vz,0,0,time,3);

  }
  
}

void SDM::SaveBoundaries_S(int itime){

  int time;
  
      time = itime - nsteps * int(itime/nsteps);

  if (Nsdm.x == 0){

    boundX(sxx,0,0,time,4);
    boundX(syy,0,0,time,5);
    boundX(szz,0,0,time,6);
    boundX(sxy,0,0,time,7);
    boundX(sxz,0,0,time,8);
    boundX(syz,0,0,time,9);

  }else if (Nsdm.x == NumSubDom.x-1) {

    boundX(sxx,1,0,time,4);
    boundX(syy,1,0,time,5);
    boundX(szz,1,0,time,6);
    boundX(sxy,1,0,time,7);
    boundX(sxz,1,0,time,8);
    boundX(syz,1,0,time,9);
  }

    if (Nsdm.y == 0){

    boundY(sxx,0,0,time,4);
    boundY(syy,0,0,time,5);
    boundY(szz,0,0,time,6);
    boundY(sxy,0,0,time,7);
    boundY(sxz,0,0,time,8);
    boundY(syz,0,0,time,9);

  }else if (Nsdm.y == NumSubDom.y-1) {

    boundY(sxx,1,0,time,4);
    boundY(syy,1,0,time,5);
    boundY(szz,1,0,time,6);
    boundY(sxy,1,0,time,7);
    boundY(sxz,1,0,time,8);
    boundY(syz,1,0,time,9);
  }

    if (Nsdm.z == 0){

    boundZ(sxx,0,0,time,4);
    boundZ(syy,0,0,time,5);
    boundZ(szz,0,0,time,6);
    boundZ(sxy,0,0,time,7);
    boundZ(sxz,0,0,time,8);
    boundZ(syz,0,0,time,9);
  }
  
}


void SDM::LoadBoundaries_V(int itime){

   int time;
 
   time = itime - nsteps * int(itime/nsteps);
    
  
  if (Nsdm.x == 0){

    boundX(vx,0,1,time,1);
    boundX(vy,0,1,time,2);
    boundX(vz,0,1,time,3);

  }else if (Nsdm.x == NumSubDom.x-1) {

    boundX(vx,1,1,time,1);
    boundX(vy,1,1,time,2);
    boundX(vz,1,1,time,3);
  }

   if (Nsdm.y == 0){

    boundY(vx,0,1,time,1);
    boundY(vy,0,1,time,2);
    boundY(vz,0,1,time,3);

  }else if (Nsdm.y == NumSubDom.y-1) {

    boundY(vx,1,1,time,1);
    boundY(vy,1,1,time,2);
    boundY(vz,1,1,time,3);
  }
   
    if (Nsdm.z == 0){

    boundZ(vx,0,1,time,1);
    boundZ(vy,0,1,time,2);
    boundZ(vz,0,1,time,3);

  }
  
  
}

void SDM::LoadBoundaries_S(int itime){

  int time;
  
      time = itime - nsteps * int(itime/nsteps);

      //printf("%d\n",time);

  if (Nsdm.x == 0){

    boundX(sxx,0,1,time,4);
    boundX(syy,0,1,time,5);
    boundX(szz,0,1,time,6);
    boundX(sxy,0,1,time,7);
    boundX(sxz,0,1,time,8);
    boundX(syz,0,1,time,9);

  }else if (Nsdm.x == NumSubDom.x-1) {

    boundX(sxx,1,1,time,4);
    boundX(syy,1,1,time,5);
    boundX(szz,1,1,time,6);
    boundX(sxy,1,1,time,7);
    boundX(sxz,1,1,time,8);
    boundX(syz,1,1,time,9);
  }

    if (Nsdm.y == 0){

    boundY(sxx,0,1,time,4);
    boundY(syy,0,1,time,5);
    boundY(szz,0,1,time,6);
    boundY(sxy,0,1,time,7);
    boundY(sxz,0,1,time,8);
    boundY(syz,0,1,time,9);

  }else if (Nsdm.y == NumSubDom.y-1) {

    boundY(sxx,1,1,time,4);
    boundY(syy,1,1,time,5);
    boundY(szz,1,1,time,6);
    boundY(sxy,1,1,time,7);
    boundY(sxz,1,1,time,8);
    boundY(syz,1,1,time,9);
  }

    if (Nsdm.z == 0){

    boundZ(sxx,0,1,time,4);
    boundZ(syy,0,1,time,5);
    boundZ(szz,0,1,time,6);
    boundZ(sxy,0,1,time,7);
    boundZ(sxz,0,1,time,8);
    boundZ(syz,0,1,time,9);
    
  }
  
}



void SDM::WriteBoundaries(int itime,int nt){

   int nt_idx;
   int Nsubs;
   int size;
   int idx_sub;
   int offset;
   char filename[200];
  
   if (((itime % nsteps) == 0) or (itime == nt )) {

     if (itime == nt) {
       nt_idx = int(floor((itime - 1)/nsteps)) + 1 ;
     } else {
       nt_idx = (itime/nsteps);
     }

     //printf("%d\n",nt_idx);

     // LEFT
      if (Nsdm.x == 0){
    

	Nsubs = NumSubDom.y * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeY() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.y + Nsdm.z * NumSubDom.y;
	offset = nt_idx * Nsubs * size + idx_sub * size ;


	if (nt_idx == 1){
	  sprintf(filename,"temp/BLX-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&BLX);
	}


	MPI_File_write(BLX,bn_lx,size,MY_MPI_Dfloat,&status);


	if (itime == nt){
	  MPI_File_close(&BLX);
	}




	// RIGHT
      }else if (Nsdm.x == NumSubDom.x-1) {

	Nsubs = NumSubDom.y * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeY() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.y + Nsdm.z * NumSubDom.y;
	offset = nt_idx * Nsubs * size + idx_sub * size;


	if (nt_idx == 1){
	  sprintf(filename,"temp/BRX-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&BRX);
	}


	MPI_File_write(BRX,bn_rx,size,MY_MPI_Dfloat,&status);


	if (itime == nt){
	  MPI_File_close(&BRX);
	}
	

      }


      // LEFT Y
      if (Nsdm.y == 0){

	Nsubs = NumSubDom.x * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.z * NumSubDom.x;
	offset = nt_idx * Nsubs * size + idx_sub * size;

	if (nt_idx == 1){
	  sprintf(filename,"temp/BLY-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&BLY);
	}


	MPI_File_write(BLY,bn_ly,size,MY_MPI_Dfloat,&status);


	if (itime == nt){
	  MPI_File_close(&BLY);
	}

    
	// RIGHT Y
      }else if (Nsdm.y == NumSubDom.y-1) {

	Nsubs = NumSubDom.x * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.z * NumSubDom.x;
	offset = nt_idx * Nsubs * size + idx_sub * size ;



	if (nt_idx == 1){
	  sprintf(filename,"temp/BRY-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&BRY);
	}


	MPI_File_write(BRY,bn_ry,size,MY_MPI_Dfloat,&status);


	if (itime == nt){
	  MPI_File_close(&BRY);
	}
    
      }
   


      // LEFT Z
      if (Nsdm.z == 0){

	Nsubs = NumSubDom.x * NumSubDom.y; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeY() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.y * NumSubDom.x;
	offset = nt_idx * Nsubs * size + idx_sub * size ;


	if (nt_idx == 1){
	  sprintf(filename,"temp/BLZ-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&BLZ);
	}


	MPI_File_write(BLZ,bn_lz,size,MY_MPI_Dfloat,&status);


	if (itime == nt){
	  MPI_File_close(&BLZ);
	}

      }


   }
  
}



void SDM::ReadBoundaries(int itime,int nt){

   int nt_idx;
   int Nsubs;
   int size;
   int idx_sub;
   MPI_Offset offset;
   char filename[200];
  
   if (((itime % nsteps) == 0) or (itime == nt)) {

     if (itime == nt) {
       nt_idx = (int(floor((itime - 1)/nsteps)) + 1) - 1;
     } else {
       nt_idx = (itime/nsteps) - 1;
     }

     //printf("%d\n",nt_idx);


     // LEFT
      if (Nsdm.x == 0){
    

	Nsubs = NumSubDom.y * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeY() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.y + Nsdm.z * NumSubDom.y;
	offset = nt_idx * size * sizeof(Dfloat);

	//printf("%d\n",offset);


	if (itime == nt){
	  sprintf(filename,"temp/BLX-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&BLX);
	}


	MPI_File_read_at(BLX,offset,bn_lx,size,MY_MPI_Dfloat,&status);


	if (nt_idx == 0){
	  MPI_File_close(&BLX);
	}




	// RIGHT
      }else if (Nsdm.x == NumSubDom.x-1) {

	Nsubs = NumSubDom.y * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeY() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.y + Nsdm.z * NumSubDom.y;
	offset = nt_idx * size* sizeof(Dfloat);


	if (itime == nt ){
	  sprintf(filename,"temp/BRX-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&BRX);
	}


	MPI_File_read_at(BRX,offset,bn_rx,size,MY_MPI_Dfloat,&status);


	if (nt_idx == 0){
	  MPI_File_close(&BRX);
	}
	

      }


      // LEFT Y
      if (Nsdm.y == 0){

	Nsubs = NumSubDom.x * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.z * NumSubDom.x;
	offset = nt_idx * size * sizeof(Dfloat);

	if (itime == nt){
	  sprintf(filename,"temp/BLY-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&BLY);
	}


	MPI_File_read_at(BLY,offset,bn_ly,size,MY_MPI_Dfloat,&status);


	if (nt_idx == 0){
	  MPI_File_close(&BLY);
	}

    
	// RIGHT Y
      }else if (Nsdm.y == NumSubDom.y-1) {

	Nsubs = NumSubDom.x * NumSubDom.z; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeZ() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.z * NumSubDom.x;
	offset = nt_idx * size * sizeof(Dfloat);



	if (itime == nt ){
	  sprintf(filename,"temp/BRY-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&BRY);
	}


	MPI_File_read_at(BRY,offset,bn_ry,size,MY_MPI_Dfloat,&status);


	if (nt_idx == 0){
	  MPI_File_close(&BRY);
	}
    
      }
   


      // LEFT Z
      if (Nsdm.z == 0){

	Nsubs = NumSubDom.x * NumSubDom.y; 
	size = KHALO * SDMGeom->L_NodeX() * SDMGeom->L_NodeY() * 9 * nsteps;
	idx_sub = Nsdm.x + Nsdm.y * NumSubDom.x;
	offset = nt_idx * size * sizeof(Dfloat);


	if (itime == nt){
	  sprintf(filename,"temp/BLZ-%d.temp",idx_sub);
	  MPI_File_open(MPI_COMM_SELF,filename,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&BLZ);
	}


	MPI_File_read_at(BLZ,offset,bn_lz,size,MY_MPI_Dfloat,&status);


	if (nt_idx == 0){
	  MPI_File_close(&BLZ);
	}

      }


   }
  
}

void SDM::PrintInf(){
  // SOURCE INFORMATION
  sourceM->PrintInf();

  // STATION INFORMATION
  station->PrintInf();
}













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


#include "kernels.hpp"

KERNEL::KERNEL(SDM *inFWD,SDM *inADJ){
  FWD = inFWD;
  ADJ = inADJ;


  KRHO = new Dfloat [FWD->SDMGeom->HALO_Node()];
  KMU = new Dfloat [FWD->SDMGeom->HALO_Node()];
  KLAMBDA = new Dfloat [FWD->SDMGeom->HALO_Node()];
  KDEN = new Dfloat [FWD->SDMGeom->HALO_Node()];
  KVP = new Dfloat [FWD->SDMGeom->HALO_Node()];
  KVS = new Dfloat [FWD->SDMGeom->HALO_Node()];
 

  ux_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  ux_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  ux_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];

  ux_dx_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  ux_dy_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  ux_dz_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uy_dx_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uy_dy_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uy_dz_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uz_dx_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uz_dy_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];
  uz_dz_ad = new Dfloat [ADJ->SDMGeom->HALO_Node()];

  InitVar(ZERO);

}


KERNEL::~KERNEL(){
  
  delete [] KRHO;
  delete [] KMU;
  delete [] KLAMBDA;
  delete [] KDEN;
  delete [] KVP;
  delete [] KVS;


  delete [] ux_dx;
  delete [] ux_dy;
  delete [] ux_dz;
  delete [] uy_dx;
  delete [] uy_dy;
  delete [] uy_dz;
  delete [] uz_dx;
  delete [] uz_dy;
  delete [] uz_dz;

  delete [] ux_dx_ad;
  delete [] ux_dy_ad;
  delete [] ux_dz_ad;
  delete [] uy_dx_ad;
  delete [] uy_dy_ad;
  delete [] uy_dz_ad;
  delete [] uz_dx_ad;
  delete [] uz_dy_ad;
  delete [] uz_dz_ad;

  FWD = NULL;
  ADJ = NULL;
  
}


void KERNEL::InitVar(Dfloat f){

  for (int i=0; i<FWD->SDMGeom->HALO_Node(); ++i){


    KRHO[i] = f;
    KMU[i] = f;
    KLAMBDA[i] = f;
    KDEN[i] = f;
    KVP[i] = f;
    KVS[i] = f;
  

  ux_dx[i] = f;
  ux_dy[i] = f;
  ux_dz[i] = f;
  uy_dx[i] = f;
  uy_dy[i] = f;
  uy_dz[i] = f;
  uz_dx[i] = f;
  uz_dy[i] = f;
  uz_dz[i] = f;

  ux_dx_ad[i] = f;
  ux_dy_ad[i] = f;
  ux_dz_ad[i] = f;
  uy_dx_ad[i] = f;
  uy_dy_ad[i] = f;
  uy_dz_ad[i] = f;
  uz_dx_ad[i] = f;
  uz_dy_ad[i] = f;
  uz_dz_ad[i] = f;

  }
  
}


void KERNEL::Dev(Dfloat *in_vx,Dfloat *in_vy,Dfloat *in_vz,Dfloat *outx_dx,Dfloat *outx_dy,Dfloat *outx_dz, \
		  Dfloat *outy_dx,Dfloat *outy_dy,Dfloat *outy_dz,\
		  Dfloat *outz_dx,Dfloat *outz_dy,Dfloat *outz_dz){
  int az;

    if (FWD->Nsdm.z == FWD->NumSubDom.z -1){
    az = FWD->SNodeZ()-KHALO -1;
  }else{
    az = FWD->SNodeZ()-KHALO;
  }
  
  for (int k=KHALO;k<az;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){


	outx_dx[FWD->IJK(i,j,k)] = ( C1 * (in_vx[FWD->IJK(i+1,j,k)] - in_vx[FWD->IJK(i,j,k)]) - \
				   C0 * (in_vx[FWD->IJK(i+2,j,k)] - in_vx[FWD->IJK(i-1,j,k)])) / FWD->SDMGeom->Dx();

	outx_dy[FWD->IJK(i,j,k)] = ( C1 * (in_vx[FWD->IJK(i,j+1,k)] - in_vx[FWD->IJK(i,j,k)]) - \
				   C0 * (in_vx[FWD->IJK(i,j+2,k)] - in_vx[FWD->IJK(i,j-1,k)])) / FWD->SDMGeom->Dy();
	
	outx_dz[FWD->IJK(i,j,k)] = ( C1 * (in_vx[FWD->IJK(i,j,k+1)] - in_vx[FWD->IJK(i,j,k)]) - \
				   C0 * (in_vx[FWD->IJK(i,j,k+2)] - in_vx[FWD->IJK(i,j,k-1)])) / FWD->SDMGeom->Dz();
	
	
	outy_dx[FWD->IJK(i,j,k)] = ( C1 * (in_vy[FWD->IJK(i,j,k)] - in_vy[FWD->IJK(i-1,j,k)]) - \
				   C0 * (in_vy[FWD->IJK(i+1,j,k)] - in_vy[FWD->IJK(i-2,j,k)])) / FWD->SDMGeom->Dx();

	outy_dy[FWD->IJK(i,j,k)] = ( C1 * (in_vy[FWD->IJK(i,j,k)] - in_vy[FWD->IJK(i,j-1,k)]) - \
				   C0 * (in_vy[FWD->IJK(i,j+1,k)] - in_vy[FWD->IJK(i,j-2,k)])) / FWD->SDMGeom->Dy();

	outy_dz[FWD->IJK(i,j,k)] = ( C1 * (in_vy[FWD->IJK(i,j,k+1)] - in_vy[FWD->IJK(i,j,k)]) - \
				   C0 * (in_vy[FWD->IJK(i,j,k+2)] - in_vy[FWD->IJK(i,j,k-1)])) / FWD->SDMGeom->Dz();


	outz_dx[FWD->IJK(i,j,k)] = ( C1 * (in_vz[FWD->IJK(i,j,k)] - in_vz[FWD->IJK(i-1,j,k)]) - \
				   C0 * (in_vz[FWD->IJK(i+1,j,k)] - in_vz[FWD->IJK(i-2,j,k)])) / FWD->SDMGeom->Dx();

	outz_dy[FWD->IJK(i,j,k)] = ( C1 * (in_vz[FWD->IJK(i,j+1,k)] - in_vz[FWD->IJK(i,j,k)]) - \
				   C0 * (in_vz[FWD->IJK(i,j+2,k)] - in_vz[FWD->IJK(i,j-1,k)])) / FWD->SDMGeom->Dy();

	outz_dz[FWD->IJK(i,j,k)] = ( C1 * (in_vz[FWD->IJK(i,j,k)] - in_vz[FWD->IJK(i,j,k-1)]) - \
				   C0 * (in_vz[FWD->IJK(i,j,k+1)] - in_vz[FWD->IJK(i,j,k-2)])) / FWD->SDMGeom->Dz();
 

      }
    }
  }

  
  
}






void KERNEL::RHO(){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	KRHO[FWD->IJK(i,j,k)] = KRHO[FWD->IJK(i,j,k)] + (	\
	  FWD->vx[FWD->IJK(i,j,k)] * ADJ->vx[ADJ->IJK(i,j,k)] + \
	  FWD->vy[FWD->IJK(i,j,k)] * ADJ->vy[ADJ->IJK(i,j,k)] + \
	  FWD->vz[FWD->IJK(i,j,k)] * ADJ->vz[ADJ->IJK(i,j,k)] ) * FWD->dt;

      }
    }
  }

}


void KERNEL::LAMBDA(){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	KLAMBDA[FWD->IJK(i,j,k)] = KLAMBDA[FWD->IJK(i,j,k)] - ((	\
	  ux_dx[FWD->IJK(i,j,k)] + uy_dy[FWD->IJK(i,j,k)] + \
	  uz_dz[FWD->IJK(i,j,k)]) * (ux_dx_ad[ADJ->IJK(i,j,k)] +	\
				     uy_dy_ad[ADJ->IJK(i,j,k)] + uz_dz_ad[ADJ->IJK(i,j,k)])) * FWD->dt;
					    
	  
      }
    }
  }

}


void KERNEL::MU(){
  Dfloat e[3][3],e_ad[3][3],buff;
  Dfloat isot;


  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){
  


	e[0][0]= ux_dx[FWD->IJK(i,j,k)];
	e[0][1]= 0.5 * (ux_dy[FWD->IJK(i,j,k)] + uy_dx[FWD->IJK(i,j,k)]);
	e[0][2]= 0.5 * (ux_dz[FWD->IJK(i,j,k)] + uz_dx[FWD->IJK(i,j,k)]);
	e[1][0]= e[0][1];
	e[1][1]= uy_dy[FWD->IJK(i,j,k)];
	e[1][2]= 0.5 * (uy_dz[FWD->IJK(i,j,k)] + uz_dy[FWD->IJK(i,j,k)]);
	e[2][0]= e[0][2];
	e[2][1]= e[1][2];
	e[2][2]= uz_dz[FWD->IJK(i,j,k)];

	isot = (e[0][0] + e[1][1] + e[2][2]) / 3.0;


	e[0][0] -= isot;
	e[1][1] -= isot;
	e[2][2] -= isot;
	
	e_ad[0][0]= ux_dx_ad[ADJ->IJK(i,j,k)];
	e_ad[0][1]= 0.5 * (ux_dy_ad[ADJ->IJK(i,j,k)] + uy_dx_ad[ADJ->IJK(i,j,k)]);
	e_ad[0][2]= 0.5 * (ux_dz_ad[ADJ->IJK(i,j,k)] + uz_dx_ad[ADJ->IJK(i,j,k)]);
	e_ad[1][0]= e_ad[0][1];
	e_ad[1][1]= uy_dy_ad[ADJ->IJK(i,j,k)];
	e_ad[1][2]= 0.5 * (uy_dz_ad[ADJ->IJK(i,j,k)] + uz_dy_ad[ADJ->IJK(i,j,k)]);
	e_ad[2][0]= e_ad[0][2];
	e_ad[2][1]= e_ad[1][2];
	e_ad[2][2]= uz_dz_ad[ADJ->IJK(i,j,k)];

	isot = (e_ad[0][0] + e_ad[1][1] + e_ad[2][2]) / 3.0;


	e_ad[0][0] -= isot;
	e_ad[1][1] -= isot;
	e_ad[2][2] -= isot;


	buff = 0.0;
	
	for (int jj = 0;jj<3;++jj){
	  for (int ii = 0;ii<3;++ii){

	    buff += e[jj][ii] * e_ad[jj][ii];

	  }
	}


	KMU[FWD->IJK(i,j,k)] = KMU[FWD->IJK(i,j,k)] - 2.0 * FWD->dt * buff;
	


      }
    }
  }

}

void KERNEL::CALC(){




  Dev(FWD->ux,FWD->uy,FWD->uz,ux_dx,ux_dy,ux_dz,	\
      uy_dx,uy_dy,uy_dz,\
      uz_dx,uz_dy,uz_dz);

  Dev(ADJ->ux,ADJ->uy,ADJ->uz,ux_dx_ad,ux_dy_ad,ux_dz_ad,	\
      uy_dx_ad,uy_dy_ad,uy_dz_ad,\
      uz_dx_ad,uz_dy_ad,uz_dz_ad);

  
  
  RHO();
  MU();
  LAMBDA();
}


void  KERNEL::KERNELS(){
   Dfloat kappa,mu,rho;
   Dfloat parmu,parkap;
   /*
Kernel Parametrization (Tromp,2005)
k bulk modulus
rho density
vp P velocity
vs S velocity
    */
 
    
    for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	  rho = FWD->rho[FWD->IJK(i,j,k)] * KRHO[FWD->IJK(i,j,k)];

	  parkap = (FWD->lamb[FWD->IJK(i,j,k)]		\
		    + (2.0/3.0) * FWD->mu[FWD->IJK(i,j,k)]);

	  parmu = FWD->mu[FWD->IJK(i,j,k)];
	  
	  
	  kappa = parkap *  KLAMBDA[FWD->IJK(i,j,k)];

	  mu = parmu * KMU[FWD->IJK(i,j,k)];

	  

	  
	  KDEN[FWD->IJK(i,j,k)] = rho + kappa + mu;

	  KVS[FWD->IJK(i,j,k)] = 2.0 * (mu - (4.0 / 3.0) * (parmu / parkap) * kappa);

	  KVP[FWD->IJK(i,j,k)] = 2.0 * ( 1.0 + (4.0 / 3.0) * (parmu / parkap)) * kappa ;
	    
	  
	}
      }
    }

}



void KERNEL::GET_K(Dfloat *KR,Dfloat *KP,Dfloat *KS){


  for (int i=0;i<FWD->SDMGeom->HALO_Node();i++){

     KR[i] = KDEN[i];
     KP[i] = KVP[i];
     KS[i] = KVS[i];

  }
}

  
  

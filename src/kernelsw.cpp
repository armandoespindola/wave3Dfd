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


#include "kernelsw.hpp"

KERNELW::KERNELW(DFT *inFWD,DFT *inADJ,SDM *sdm){

  fwdDF = inFWD;
  adjDF = inADJ;
  FWD = sdm;


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

  ux_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  ux_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  ux_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uy_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  uz_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];

  iKRHO = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iKMU = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iKLAMBDA = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iKDEN = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iKVP = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iKVS = new Dfloat [FWD->SDMGeom->HALO_Node()];
 

  iux_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iux_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iux_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dx = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dy = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dz = new Dfloat [FWD->SDMGeom->HALO_Node()];

  iux_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iux_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iux_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuy_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dx_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dy_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];
  iuz_dz_ad = new Dfloat [FWD->SDMGeom->HALO_Node()];

  iPcondA = new Dfloat [FWD->SDMGeom->HALO_Node()];
  PcondA = new Dfloat [FWD->SDMGeom->HALO_Node()];

  iPcondB = new Dfloat [FWD->SDMGeom->HALO_Node()];
  PcondB = new Dfloat [FWD->SDMGeom->HALO_Node()];

  InitVar(ZERO);

}


KERNELW::~KERNELW(){
  
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

  delete [] iKRHO;
  delete [] iKMU;
  delete [] iKLAMBDA;
  delete [] iKDEN;
  delete [] iKVP;
  delete [] iKVS;


  delete [] iux_dz;
  delete [] iux_dy;
  delete [] iux_dx;
  delete [] iuy_dx;
  delete [] iuy_dy;
  delete [] iuy_dz;
  delete [] iuz_dx;
  delete [] iuz_dy;
  delete [] iuz_dz;

  delete [] iux_dx_ad;
  delete [] iux_dy_ad;
  delete [] iux_dz_ad;
  delete [] iuy_dx_ad;
  delete [] iuy_dy_ad;
  delete [] iuy_dz_ad;
  delete [] iuz_dx_ad;
  delete [] iuz_dy_ad;
  delete [] iuz_dz_ad;

  delete [] iPcondA;
  delete [] PcondA;

  delete [] iPcondB;
  delete [] PcondB;

  fwdDF = NULL;
  adjDF = NULL;
  FWD = NULL;
  
}


void KERNELW::InitVar(Dfloat f){

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

   iKRHO[i] = f;
    iKMU[i] = f;
    iKLAMBDA[i] = f;
    iKDEN[i] = f;
    iKVP[i] = f;
    iKVS[i] = f;
  

  iux_dx[i] = f;
  iux_dy[i] = f;
  iux_dz[i] = f;
  iuy_dx[i] = f;
  iuy_dy[i] = f;
  iuy_dz[i] = f;
  iuz_dx[i] = f;
  iuz_dy[i] = f;
  iuz_dz[i] = f;

  iux_dx_ad[i] = f;
  iux_dy_ad[i] = f;
  iux_dz_ad[i] = f;
  iuy_dx_ad[i] = f;
  iuy_dy_ad[i] = f;
  iuy_dz_ad[i] = f;
  iuz_dx_ad[i] = f;
  iuz_dy_ad[i] = f;
  iuz_dz_ad[i] = f;

  iPcondA[i] = f;
  PcondA[i] = f;

  iPcondB[i] = f;
  PcondB[i] = f;

  }
  
}
void KERNELW::DevX(Dfloat **in_var,int l,Dfloat *out_var){
int az;

  if (FWD->Nsdm.z == FWD->NumSubDom.z -1){
    az = FWD->SNodeZ()-KHALO -1;
  }else{
    az = FWD->SNodeZ()-KHALO;
  }
  
  for (int k=KHALO;k<az;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){


	//	if ((i >= 2*KHALO) && (i < FWD->SNodeX()-2*KHALO)){
	  
	  out_var[FWD->IJK(i,j,k)] = ( (1.0 / 12.0) * in_var[l][FWD->IJK(i-2,j,k)] - \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i-1,j,k)] + \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i+1,j,k)] - \
				(1.0 / 12.0) * in_var[l][FWD->IJK(i+2,j,k)] ) / FWD->SDMGeom->Dx();
	  //	}


	//RIGHT SIDE

	/*	if (i < 2*KHALO) {

	  out_var[FWD->IJK(i,j,k)] = ( - (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] + \
				 4.0 * in_var[l][FWD->IJK(i+1,j,k)] - \
				 3.0 * in_var[l][FWD->IJK(i+2,j,k)] + \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i+3,j,k)] -
				  (1.0 / 4.0) * in_var[l][FWD->IJK(i+4,j,k)] ) / FWD->SDMGeom->Dx();
	}
	*/
	// LEFT SIDE
	/*	
	if (i >= FWD->SNodeX()-2*KHALO){

	  out_var[FWD->IJK(i,j,k)] = (  (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] - \
				 4.0 * in_var[l][FWD->IJK(i-1,j,k)] + \
				 3.0 * in_var[l][FWD->IJK(i-2,j,k)] - \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i-3,j,k)] +
				   (1.0 / 4.0) * in_var[l][FWD->IJK(i-4,j,k)] ) / FWD->SDMGeom->Dx();

	}
	
	*/
	}
      }
  }

  

  
}


void KERNELW::DevY(Dfloat **in_var,int l,Dfloat *out_var){
  int az;

  if (FWD->Nsdm.z == FWD->NumSubDom.z -1){
    az = FWD->SNodeZ()-KHALO -1;
  }else{
    az = FWD->SNodeZ()-KHALO;
  }
    for (int k=KHALO;k<az;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){


	//	if ((j >= 2*KHALO) && (j < FWD->SNodeY()-2*KHALO)){
	  
	out_var[FWD->IJK(i,j,k)] = ( (1.0 / 12.0) * in_var[l][FWD->IJK(i,j-2,k)] - \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i,j-1,k)] + \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i,j+1,k)] - \
				(1.0 / 12.0) * in_var[l][FWD->IJK(i,j+2,k)] ) / FWD->SDMGeom->Dy();
	//	}


	//RIGHT SIDE

	/*if (j < 2*KHALO) {

	  out_var[FWD->IJK(i,j,k)] = ( - (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] + \
				 4.0 * in_var[l][FWD->IJK(i,j+1,k)] - \
				 3.0 * in_var[l][FWD->IJK(i,j+2,k)] + \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i,j+3,k)] -
				  (1.0 / 4.0) * in_var[l][FWD->IJK(i,j+4,k)] ) / FWD->SDMGeom->Dy();
	}
	*/
	// LEFT SIDE
	/*
	if (j >= FWD->SNodeY()-2*KHALO){

	  out_var[FWD->IJK(i,j,k)] = (  (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] - \
				 4.0 * in_var[l][FWD->IJK(i,j-1,k)] + \
				 3.0 * in_var[l][FWD->IJK(i,j-2,k)] - \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i,j-3,k)] +
				  (1.0 / 4.0) * in_var[l][FWD->IJK(i,j-4,k)] ) / FWD->SDMGeom->Dy();

	}
	*/

	}
      }
    }

    
}



void KERNELW::DevZ(Dfloat **in_var,int l,Dfloat *out_var){
  int az;

  if (FWD->Nsdm.z == FWD->NumSubDom.z -1){
    az = FWD->SNodeZ()-KHALO -1;
  }else{
    az = FWD->SNodeZ()-KHALO;
  }
    

  for (int k=KHALO;k<az;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){


	//	if ((k >= 2*KHALO) && (k < FWD->SNodeZ()-2*KHALO-1)){
	  
	out_var[FWD->IJK(i,j,k)] = ( (1.0 / 12.0) * in_var[l][FWD->IJK(i,j,k-2)] - \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i,j,k-1)] + \
				(2.0 / 3.0) * in_var[l][FWD->IJK(i,j,k+1)] - \
				(1.0 / 12.0) * in_var[l][FWD->IJK(i,j,k+2)] ) / FWD->SDMGeom->Dz();
	//}


	//RIGHT SIDE

	/*if (k < 2*KHALO) {

	  out_var[FWD->IJK(i,j,k)] = ( - (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] + \
				 4.0 * in_var[l][FWD->IJK(i,j,k+1)] - \
				 3.0 * in_var[l][FWD->IJK(i,j,k+2)] + \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i,j,k+3)] -
				  (1.0 / 4.0) * in_var[l][FWD->IJK(i,j,k+4)] ) / FWD->SDMGeom->Dz();
	}
	*/
	// LEFT SIDE
	/*
	if (k >= FWD->SNodeZ()-2*KHALO-1){

	  out_var[FWD->IJK(i,j,k)] = (  (25.0 / 12.0) * in_var[l][FWD->IJK(i,j,k)] - \
				 4.0 * in_var[l][FWD->IJK(i,j,k-1)] + \
				 3.0 * in_var[l][FWD->IJK(i,j,k-2)] - \
				(4.0 / 3.0) * in_var[l][FWD->IJK(i,j,k-3)] +
				  (1.0 / 4.0) * in_var[l][FWD->IJK(i,j,k-4)] ) / FWD->SDMGeom->Dz();

	}
	*/

	}
      }
    }

}



void KERNELW::RHO(int l){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	KRHO[FWD->IJK(i,j,k)] = 1.0  * \
	  ((fwdDF->Fvx[l][FWD->IJK(i,j,k)] * adjDF->Fvx[l][FWD->IJK(i,j,k)]   + \
	   fwdDF->Fvy[l][FWD->IJK(i,j,k)] * adjDF->Fvy[l][FWD->IJK(i,j,k)]   + \
	    fwdDF->Fvz[l][FWD->IJK(i,j,k)] * adjDF->Fvz[l][FWD->IJK(i,j,k)] )   + \
	   (fwdDF->iFvx[l][FWD->IJK(i,j,k)] * adjDF->iFvx[l][FWD->IJK(i,j,k)]   + \
	   fwdDF->iFvy[l][FWD->IJK(i,j,k)] * adjDF->iFvy[l][FWD->IJK(i,j,k)]   + \
	    fwdDF->iFvz[l][FWD->IJK(i,j,k)] * adjDF->iFvz[l][FWD->IJK(i,j,k)]));

	PcondB[FWD->IJK(i,j,k)] += KRHO[FWD->IJK(i,j,k)] * fwdDF->freq[l] * fwdDF->freq[l] * 4.0 * pi * pi  ;

  PcondA[FWD->IJK(i,j,k)] +=  -4.0 * pi * pi * fwdDF->freq[l] * fwdDF->freq[l] * \
    ((fwdDF->Fvx[l][FWD->IJK(i,j,k)] * fwdDF->Fvx[l][FWD->IJK(i,j,k)]   + \
     fwdDF->Fvy[l][FWD->IJK(i,j,k)] * fwdDF->Fvy[l][FWD->IJK(i,j,k)]   + \
      fwdDF->Fvz[l][FWD->IJK(i,j,k)] * fwdDF->Fvz[l][FWD->IJK(i,j,k)] )   - \
     (fwdDF->iFvx[l][FWD->IJK(i,j,k)] * fwdDF->iFvx[l][FWD->IJK(i,j,k)]   + \
     fwdDF->iFvy[l][FWD->IJK(i,j,k)] * fwdDF->iFvy[l][FWD->IJK(i,j,k)]   + \
      fwdDF->iFvz[l][FWD->IJK(i,j,k)] * fwdDF->iFvz[l][FWD->IJK(i,j,k)])); 

      }
    }
  }

}

/*
void KERNELW::iRHO(int l){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	iKRHO[FWD->IJK(i,j,k)] =  fwdDF->freq[l] * fwdDF->freq[l] * \
	  ((fwdDF->iFux[l][FWD->IJK(i,j,k)] * adjDF->Fux[l][FWD->IJK(i,j,k)]   + \
	   fwdDF->iFuy[l][FWD->IJK(i,j,k)] * adjDF->Fuy[l][FWD->IJK(i,j,k)]   + \
	    fwdDF->iFuz[l][FWD->IJK(i,j,k)] * adjDF->Fuz[l][FWD->IJK(i,j,k)] )   - \
	   (fwdDF->Fux[l][FWD->IJK(i,j,k)] * adjDF->iFux[l][FWD->IJK(i,j,k)]   + \
	   fwdDF->Fuy[l][FWD->IJK(i,j,k)] * adjDF->iFuy[l][FWD->IJK(i,j,k)]   + \
	    fwdDF->Fuz[l][FWD->IJK(i,j,k)] * adjDF->iFuz[l][FWD->IJK(i,j,k)])) * FWD->dt;

	//iPcondB[FWD->IJK(i,j,k)] += iKRHO[FWD->IJK(i,j,k)] * fwdDF->freq[l] * fwdDF->freq[l];

      }
    }
  }

}
*/

void KERNELW::LAMBDA(int l){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	KLAMBDA[FWD->IJK(i,j,k)] = ((   \
		(ux_dx[FWD->IJK(i,j,k)] + uy_dy[FWD->IJK(i,j,k)] + \
	  uz_dz[FWD->IJK(i,j,k)] ) * (ux_dx_ad[FWD->IJK(i,j,k)] + \
				      uy_dy_ad[FWD->IJK(i,j,k)] + uz_dz_ad[FWD->IJK(i,j,k)])) + \
	    ((iux_dx[FWD->IJK(i,j,k)] + iuy_dy[FWD->IJK(i,j,k)] + \
	  iuz_dz[FWD->IJK(i,j,k)] ) * (iux_dx_ad[FWD->IJK(i,j,k)] + \
				       iuy_dy_ad[FWD->IJK(i,j,k)] + iuz_dz_ad[FWD->IJK(i,j,k)]))) \
	  * - (1.0 / (fwdDF->freq[l] * fwdDF->freq[l] * 4.0 * pi * pi)) ;
					    
	  
      }
    }
  }

}

/*
void KERNELW::iLAMBDA(){

  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	iKLAMBDA[FWD->IJK(i,j,k)] = ((   \
		(iux_dx[FWD->IJK(i,j,k)] + iuy_dy[FWD->IJK(i,j,k)] + \
	  iuz_dz[FWD->IJK(i,j,k)] ) * (ux_dx_ad[FWD->IJK(i,j,k)] + \
	    uy_dy_ad[FWD->IJK(i,j,k)] + uz_dz_ad[FWD->IJK(i,j,k)])) - \
	    ((ux_dx[FWD->IJK(i,j,k)] + uy_dy[FWD->IJK(i,j,k)] + \
	  uz_dz[FWD->IJK(i,j,k)] ) * (iux_dx_ad[FWD->IJK(i,j,k)] + \
	    iuy_dy_ad[FWD->IJK(i,j,k)] + iuz_dz_ad[FWD->IJK(i,j,k)]))) * FWD->dt;
					    
	  
      }
    }
  }

}

*/

void KERNELW::MU(int l){
  Dfloat e[3][3],e_ad[3][3],buff;
  Dfloat ie[3][3],ie_ad[3][3],ibuff;
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


	e_ad[0][0]= ux_dx_ad[FWD->IJK(i,j,k)];
	e_ad[0][1]= 0.5 * (ux_dy_ad[FWD->IJK(i,j,k)] + uy_dx_ad[FWD->IJK(i,j,k)]);
	e_ad[0][2]= 0.5 * (ux_dz_ad[FWD->IJK(i,j,k)] + uz_dx_ad[FWD->IJK(i,j,k)]);
	e_ad[1][0]= e_ad[0][1];
	e_ad[1][1]= uy_dy_ad[FWD->IJK(i,j,k)];
	e_ad[1][2]= 0.5 * (uy_dz_ad[FWD->IJK(i,j,k)] + uz_dy_ad[FWD->IJK(i,j,k)]);
	e_ad[2][0]= e_ad[0][2];
	e_ad[2][1]= e_ad[1][2];
	e_ad[2][2]= uz_dz_ad[FWD->IJK(i,j,k)];

	isot = (e_ad[0][0] + e_ad[1][1] + e_ad[2][2]) / 3.0;

	e_ad[0][0] -= isot;
	e_ad[1][1] -= isot;
	e_ad[2][2] -= isot;


	buff = 0.0;

	ie[0][0]= iux_dx[FWD->IJK(i,j,k)];
	ie[0][1]= 0.5 * (iux_dy[FWD->IJK(i,j,k)] + iuy_dx[FWD->IJK(i,j,k)]);
	ie[0][2]= 0.5 * (iux_dz[FWD->IJK(i,j,k)] + iuz_dx[FWD->IJK(i,j,k)]);
	ie[1][0]= ie[0][1];
	ie[1][1]= iuy_dy[FWD->IJK(i,j,k)];
	ie[1][2]= 0.5 * (iuy_dz[FWD->IJK(i,j,k)] + iuz_dy[FWD->IJK(i,j,k)]);
	ie[2][0]= ie[0][2];
	ie[2][1]= ie[1][2];
	ie[2][2]= iuz_dz[FWD->IJK(i,j,k)];

	isot = (ie[0][0] + ie[1][1] + ie[2][2]) / 3.0;

	ie[0][0] -= isot;
	ie[1][1] -= isot;
	ie[2][2] -= isot;


	ie_ad[0][0]= iux_dx_ad[FWD->IJK(i,j,k)];
	ie_ad[0][1]= 0.5 * (iux_dy_ad[FWD->IJK(i,j,k)] + iuy_dx_ad[FWD->IJK(i,j,k)]);
	ie_ad[0][2]= 0.5 * (iux_dz_ad[FWD->IJK(i,j,k)] + iuz_dx_ad[FWD->IJK(i,j,k)]);
	ie_ad[1][0]= ie_ad[0][1];
	ie_ad[1][1]= iuy_dy_ad[FWD->IJK(i,j,k)];
	ie_ad[1][2]= 0.5 * (iuy_dz_ad[FWD->IJK(i,j,k)] + iuz_dy_ad[FWD->IJK(i,j,k)]);
	ie_ad[2][0]= ie_ad[0][2];
	ie_ad[2][1]= ie_ad[1][2];
	ie_ad[2][2]= iuz_dz_ad[FWD->IJK(i,j,k)];

	isot = (ie_ad[0][0] + ie_ad[1][1] + ie_ad[2][2]) / 3.0;

	ie_ad[0][0] -= isot;
	ie_ad[1][1] -= isot;
	ie_ad[2][2] -= isot;


	ibuff = 0.0;
	
	for (int jj = 0;jj<3;++jj){
	  for (int ii = 0;ii<3;++ii){

	    buff += e[jj][ii] * e_ad[jj][ii];
	    ibuff += ie[jj][ii] * ie_ad[jj][ii];

	  }
	}


	KMU[FWD->IJK(i,j,k)] = (-2.0) * (buff + ibuff) *		\
	   (1.0 / (fwdDF->freq[l] * fwdDF->freq[l] * 4.0 * pi * pi)) ;
	


      }
    }
  }

}

/*
void KERNELW::iMU(){
  Dfloat e[3][3],e_ad[3][3],buff;
  Dfloat ie[3][3],ie_ad[3][3],ibuff;


  for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
    for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
      for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){
  


	e[0][0]= 2.0 * ux_dx[FWD->IJK(i,j,k)];
	e[0][1]= ux_dy[FWD->IJK(i,j,k)] + uy_dx[FWD->IJK(i,j,k)];
	e[0][2]= ux_dz[FWD->IJK(i,j,k)] + uz_dx[FWD->IJK(i,j,k)];
	e[1][0]= e[0][1];
	e[1][1]= 2.0 * uy_dy[FWD->IJK(i,j,k)];
	e[1][2]= uy_dz[FWD->IJK(i,j,k)] + uz_dy[FWD->IJK(i,j,k)];
	e[2][0]= e[0][2];
	e[2][1]= e[1][2];
	e[2][2]= 2.0 * uz_dz[FWD->IJK(i,j,k)];


	e_ad[0][0]= 2.0 * ux_dx_ad[FWD->IJK(i,j,k)];
	e_ad[0][1]= ux_dy_ad[FWD->IJK(i,j,k)] + uy_dx_ad[FWD->IJK(i,j,k)];
	e_ad[0][2]= ux_dz_ad[FWD->IJK(i,j,k)] + uz_dx_ad[FWD->IJK(i,j,k)];
	e_ad[1][0]= e_ad[0][1];
	e_ad[1][1]= 2.0 * uy_dy_ad[FWD->IJK(i,j,k)];
	e_ad[1][2]= uy_dz_ad[FWD->IJK(i,j,k)] + uz_dy_ad[FWD->IJK(i,j,k)];
	e_ad[2][0]= e_ad[0][2];
	e_ad[2][1]= e_ad[1][2];
	e_ad[2][2]= 2.0 * uz_dz_ad[FWD->IJK(i,j,k)];


	buff = 0.0;

	ie[0][0]= 2.0 * iux_dx[FWD->IJK(i,j,k)];
	ie[0][1]= iux_dy[FWD->IJK(i,j,k)] + iuy_dx[FWD->IJK(i,j,k)];
	ie[0][2]= iux_dz[FWD->IJK(i,j,k)] + iuz_dx[FWD->IJK(i,j,k)];
	ie[1][0]= ie[0][1];
	ie[1][1]= 2.0 * iuy_dy[FWD->IJK(i,j,k)];
	ie[1][2]= iuy_dz[FWD->IJK(i,j,k)] + iuz_dy[FWD->IJK(i,j,k)];
	ie[2][0]= ie[0][2];
	ie[2][1]= ie[1][2];
	ie[2][2]= 2.0 * iuz_dz[FWD->IJK(i,j,k)];


	ie_ad[0][0]= 2.0 * iux_dx_ad[FWD->IJK(i,j,k)];
	ie_ad[0][1]= iux_dy_ad[FWD->IJK(i,j,k)] + iuy_dx_ad[FWD->IJK(i,j,k)];
	ie_ad[0][2]= iux_dz_ad[FWD->IJK(i,j,k)] + iuz_dx_ad[FWD->IJK(i,j,k)];
	ie_ad[1][0]= ie_ad[0][1];
	ie_ad[1][1]= 2.0 * iuy_dy_ad[FWD->IJK(i,j,k)];
	ie_ad[1][2]= iuy_dz_ad[FWD->IJK(i,j,k)] + iuz_dy_ad[FWD->IJK(i,j,k)];
	ie_ad[2][0]= ie_ad[0][2];
	ie_ad[2][1]= ie_ad[1][2];
	ie_ad[2][2]= 2.0 * iuz_dz_ad[FWD->IJK(i,j,k)];


	ibuff = 0.0;
	
	for (int jj = 0;jj<3;++jj){
	  for (int ii = 0;ii<3;++ii){

	    buff += ie[jj][ii] * e_ad[jj][ii];
	    ibuff += e[jj][ii] * ie_ad[jj][ii];

	  }
	}


	iKMU[FWD->IJK(i,j,k)] = 0.5 * FWD->dt * (buff - ibuff);
	


      }
    }
  }

}

*/

void KERNELW::CALC(int i){

  //InitVar(ZERO);
 
  DevX(fwdDF->Fvx,i,ux_dx);
  DevY(fwdDF->Fvx,i,ux_dy);
  DevZ(fwdDF->Fvx,i,ux_dz);
  
  DevX(fwdDF->Fvy,i,uy_dx);
  DevY(fwdDF->Fvy,i,uy_dy);
  DevZ(fwdDF->Fvy,i,uy_dz);
  
  DevX(fwdDF->Fvz,i,uz_dx);
  DevY(fwdDF->Fvz,i,uz_dy);
  DevZ(fwdDF->Fvz,i,uz_dz);

  DevX(adjDF->Fvx,i,ux_dx_ad);
  DevY(adjDF->Fvx,i,ux_dy_ad);
  DevZ(adjDF->Fvx,i,ux_dz_ad);
  
  DevX(adjDF->Fvy,i,uy_dx_ad);
  DevY(adjDF->Fvy,i,uy_dy_ad);
  DevZ(adjDF->Fvy,i,uy_dz_ad);
  
  DevX(adjDF->Fvz,i,uz_dx_ad);
  DevY(adjDF->Fvz,i,uz_dy_ad);
  DevZ(adjDF->Fvz,i,uz_dz_ad);


  DevX(fwdDF->iFvx,i,iux_dx);
  DevY(fwdDF->iFvx,i,iux_dy);
  DevZ(fwdDF->iFvx,i,iux_dz);
  
  DevX(fwdDF->iFvy,i,iuy_dx);
  DevY(fwdDF->iFvy,i,iuy_dy);
  DevZ(fwdDF->iFvy,i,iuy_dz);
  
  DevX(fwdDF->iFvz,i,iuz_dx);
  DevY(fwdDF->iFvz,i,iuz_dy);
  DevZ(fwdDF->iFvz,i,iuz_dz);

  DevX(adjDF->iFvx,i,iux_dx_ad);
  DevY(adjDF->iFvx,i,iux_dy_ad);
  DevZ(adjDF->iFvx,i,iux_dz_ad);
  
  DevX(adjDF->iFvy,i,iuy_dx_ad);
  DevY(adjDF->iFvy,i,iuy_dy_ad);
  DevZ(adjDF->iFvy,i,iuy_dz_ad);
  
  DevX(adjDF->iFvz,i,iuz_dx_ad);
  DevY(adjDF->iFvz,i,iuz_dy_ad);
  DevZ(adjDF->iFvz,i,iuz_dz_ad);
  
  
  RHO(i);
  MU(i);
  LAMBDA(i);
  KERNELS();
 
  //iRHO(i);
  //iMU();
  //iLAMBDA();
  //iKERNELS();
  
  
  
}


void  KERNELW::KERNELS(){
   Dfloat kappa,mu,rho;
   Dfloat parmu,parkap;
   /*
Kernel Parametrization (Tromp,2005)
k bulk modulus
rho density
vp P velocity
vs S velocity
    */

   int az;

  if (FWD->Nsdm.z == FWD->NumSubDom.z -1){
    az = FWD->SNodeZ()-KHALO -1;
  }else{
    az = FWD->SNodeZ()-KHALO;
  }
 
    
    for (int k=KHALO;k<az;k++){
      for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	  rho = FWD->rho[FWD->IJK(i,j,k)] * KRHO[FWD->IJK(i,j,k)];

	  parkap = (FWD->lamb[FWD->IJK(i,j,k)]		\
		    + (2.0/3.0) * FWD->mu[FWD->IJK(i,j,k)]);

	  parmu = FWD->mu[FWD->IJK(i,j,k)];
	  
	  
	  kappa = parkap *  KLAMBDA[FWD->IJK(i,j,k)];

	  mu = parmu * KMU[FWD->IJK(i,j,k)];

	  KDEN[FWD->IJK(i,j,k)] += rho + kappa + mu;

	  KVS[FWD->IJK(i,j,k)] += 2.0 * (mu - (4.0 / 3.0) * (parmu / parkap) * kappa);

	  KVP[FWD->IJK(i,j,k)] += 2.0 * ( 1.0 + (4.0 / 3.0) * (parmu / parkap)) * kappa ;
	    
	  
	}
      }
    }

}



void KERNELW::GET_K(Dfloat *KR,Dfloat *KP,Dfloat *KS,Dfloat *KPA,Dfloat *KPB){


  for (int i=0;i<FWD->SDMGeom->HALO_Node();i++){

     KR[i] = KDEN[i];
     KP[i] = KVP[i];
     KS[i] = KVS[i];
     KPA[i] = PcondA[i];
     KPB[i] = PcondB[i];

  }
}








/*
void  KERNELW::iKERNELS(){
   Dfloat vp,vs,rho;

 
    
    for (int k=KHALO;k<FWD->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<FWD->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<FWD->SNodeX()-KHALO;i++){

	  rho = FWD->rho[FWD->IJK(i,j,k)];
	  
	  vp = (FWD->lamb[FWD->IJK(i,j,k)] \
		+ 2.0 * FWD->mu[FWD->IJK(i,j,k)]) / rho;

	  vs = FWD->mu[FWD->IJK(i,j,k)]  / rho;

	  
	  iKDEN[FWD->IJK(i,j,k)] += iKRHO[FWD->IJK(i,j,k)] +	\
	    (vp - 2.0 * vs) * iKLAMBDA[FWD->IJK(i,j,k)] + \
	    vs * iKMU[FWD->IJK(i,j,k)];

	  iKVS[FWD->IJK(i,j,k)] += 2.0 * rho * sqrt(vs) * iKMU[FWD->IJK(i,j,k)] - \
	    4.0 * rho * sqrt(vs) * iKLAMBDA[FWD->IJK(i,j,k)];

	  iKVP[FWD->IJK(i,j,k)] += 2.0 * rho * sqrt(vp) * iKLAMBDA[FWD->IJK(i,j,k)];
	    
	  
	}
      }
    }

}

*/

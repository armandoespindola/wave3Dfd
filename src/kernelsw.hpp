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


#ifndef __KERNELW__
#define __KERNELW__

#include "definitions.hpp"
#include "sdm.hpp"
#include "fdt.hpp"

class KERNELW {

protected:
  // KRHO Kernel density
  // KMU Kernel Mu
  // KLAMBDA Kernel Lambda
  // Fitchner (2010)
 
  DFT *fwdDF,*adjDF;
  SDM *FWD;
  
  Dfloat *ux_dx,*ux_dy,*ux_dz;
  Dfloat *uy_dx,*uy_dy,*uy_dz;
  Dfloat *uz_dx,*uz_dy,*uz_dz;
  
  Dfloat *ux_dx_ad,*ux_dy_ad,*ux_dz_ad;
  Dfloat *uy_dx_ad,*uy_dy_ad,*uy_dz_ad;
  Dfloat *uz_dx_ad,*uz_dy_ad,*uz_dz_ad;


  Dfloat *iux_dx,*iux_dy,*iux_dz;
  Dfloat *iuy_dx,*iuy_dy,*iuy_dz;
  Dfloat *iuz_dx,*iuz_dy,*iuz_dz;
  
  Dfloat *iux_dx_ad,*iux_dy_ad,*iux_dz_ad;
  Dfloat *iuy_dx_ad,*iuy_dy_ad,*iuy_dz_ad;
  Dfloat *iuz_dx_ad,*iuz_dy_ad,*iuz_dz_ad;


public:


  Dfloat *KRHO,*KMU,*KLAMBDA,*KDEN,*KVP,*KVS;
  Dfloat *iKRHO,*iKMU,*iKLAMBDA,*iKDEN,*iKVP,*iKVS;

  // Preconditioner
  Dfloat *PcondA,*iPcondA;
  Dfloat *PcondB,*iPcondB;

  KERNELW(DFT *inFWD,DFT *inADJ,SDM *sdm);

  ~KERNELW();

  void InitVar(Dfloat f);

  void RHO(int i);
  void MU(int i);
  void LAMBDA(int i);
  void KERNELS();
  void CALC(int i);

  /*
  void iRHO(int i);
  void iMU();
  void iLAMBDA();
  void iKERNELS();
  */
  
  // SPATIAL DERIVATIVES (4TH ORDER)

  void Dev(int l,Dfloat **in_vx,Dfloat **in_vy,Dfloat **in_vz,Dfloat *outx_dx,Dfloat *outx_dy,Dfloat *outx_dz, \
	   Dfloat *outy_dx,Dfloat *outy_dy,Dfloat *outy_dz,\
	   Dfloat *outz_dx,Dfloat *outz_dy,Dfloat *outz_dz);
  
  void GET_K(Dfloat *KR,Dfloat *KP,Dfloat *KS,Dfloat *KPA,Dfloat *KPB);

};

#endif


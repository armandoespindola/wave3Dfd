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


#ifndef __KERNEL__
#define __KERNEL__

#include "definitions.hpp"
#include "sdm.hpp"

class KERNEL {

protected:
  // KRHO Kernel density
  // KMU Kernel Mu
  // KLAMBDA Kernel Lambda
  // Kernel Parametrization (Tromp,2005)
  // FWD Forward Propagation
  // ADJ Adjoint Propagation
  
  SDM *FWD,*ADJ;

public:

  Dfloat *ux_dx,*ux_dy,*ux_dz;
  Dfloat *uy_dx,*uy_dy,*uy_dz;
  Dfloat *uz_dx,*uz_dy,*uz_dz;
  
  Dfloat *ux_dx_ad,*ux_dy_ad,*ux_dz_ad;
  Dfloat *uy_dx_ad,*uy_dy_ad,*uy_dz_ad;
  Dfloat *uz_dx_ad,*uz_dy_ad,*uz_dz_ad;


  Dfloat *KRHO,*KMU,*KLAMBDA,*KDEN,*KVP,*KVS;

  KERNEL(SDM *inFWD,SDM *inADJ);

  ~KERNEL();

  void InitVar(Dfloat f);

  void RHO();
  void MU();
  void LAMBDA();
  void KERNELS();
  void CALC();

  // SPATIAL DERIVATIVES (4TH ORDER)

   void Dev(Dfloat *in_vx,Dfloat *in_vy,Dfloat *in_vz,Dfloat *outx_dx,Dfloat *outx_dy,Dfloat *outx_dz, \
	   Dfloat *outy_dx,Dfloat *outy_dy,Dfloat *outy_dz,\
	   Dfloat *outz_dx,Dfloat *outz_dy,Dfloat *outz_dz);


  void GET_K(Dfloat *KR,Dfloat *KP,Dfloat *KS);

};

#endif

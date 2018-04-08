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

  KERNELW(DFT *inFWD,DFT *inADJ,SDM *sdm);

  ~KERNELW();

  void InitVar(Dfloat f);

  void RHO(int i);
  void MU();
  void LAMBDA();
  void KERNELS();
  void CALC(int i);

  void iRHO(int i);
  void iMU();
  void iLAMBDA();
  void iKERNELS();

  // SPATIAL DERIVATIVES (4TH ORDER)
  
  void DevX(Dfloat *in_var,Dfloat *out_var);
  void DevY(Dfloat *in_var,Dfloat *out_var);
  void DevZ(Dfloat *in_var,Dfloat *out_var);

};

#endif


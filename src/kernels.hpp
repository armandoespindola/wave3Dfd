#ifndef __KERNEL__
#define __KERNEL__

#include "definitions.hpp"
#include "sdm.hpp"

class KERNEL {

protected:
  // KRHO Kernel density
  // KMU Kernel Mu
  // KLAMBDA Kernel Lambda
  // Fitchner (2010)
  // FWD Forward Propagation
  // ADJ Adjoint Propagation
  
  SDM *FWD,*ADJ;
  
  Dfloat *ux_dx,*ux_dy,*ux_dz;
  Dfloat *uy_dx,*uy_dy,*uy_dz;
  Dfloat *uz_dx,*uz_dy,*uz_dz;
  
  Dfloat *ux_dx_ad,*ux_dy_ad,*ux_dz_ad;
  Dfloat *uy_dx_ad,*uy_dy_ad,*uy_dz_ad;
  Dfloat *uz_dx_ad,*uz_dy_ad,*uz_dz_ad;
  int nfq;

public:


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
  
  void DevX(Dfloat *in_var,Dfloat *out_var);
  void DevY(Dfloat *in_var,Dfloat *out_var);
  void DevZ(Dfloat *in_var,Dfloat *out_var);

};

#endif

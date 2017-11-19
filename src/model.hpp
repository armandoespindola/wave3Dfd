#ifndef __MODEL__
#define __MODEL__

#include "definitions.hpp"

class MODEL {

protected:
  // GDIM GLOBAL DIMENSION NODES WITHOUT PML
  // SubDomNodeN NUMBER OF NODES SUBDOMAINS
  // HALO HALO FINITE DIFFERENCE
  // FileName NAME OF FILE 
  VecI GDim;
  VecI SubDomNodeN;
  Dfloat *Model,*ModelPML;
  VecI NDT,NDL;
  Dfloat *mu,*rho,*lamb;
  Dfloat *muL,*rhoL,*lambL;
  VecI Nsub;

public:

  MODEL(char *FileP,char *FileS,char *FileR, VecI iGDim, VecI iSubDomNodeN);

  ~MODEL();

  void PML_MODEL();

  void SubModel(VecI NumSub, Dfloat *SubR, Dfloat *SubM, Dfloat *SubL);
  
  
  

};


#endif

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

  MODEL(std::string FileP,std::string FileS,std::string FileR, VecI iGDim, VecI iSubDomNodeN);

  ~MODEL();

  void PML_MODEL();

  void SubModel(VecI NumSub, Dfloat *SubR, Dfloat *SubM, Dfloat *SubL);

// RETURN NUMBER OF TOTAL NODES WITH (HALO OR PML)
  inline int HALO_Node() { return (NDL.x * NDL.y * NDL.z); }

  // NUMBER OF NODES IN X DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeX() { return NDL.x; }

  // NUMBER OF NODES IN Y DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeY() { return NDL.y; }

  // NUMBER OF NODES IN Z DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeZ() { return NDL.z; }
  
  

};


#endif

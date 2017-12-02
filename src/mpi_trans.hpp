#ifndef __MPITRANS__
#define __MPITRANS__

#include "definitions.hpp"
#include "geometry3D.hpp"
#include "model.hpp"
#include "show.hpp"
#include "sdm.hpp"


class MPI_DATA {
protected:
	SDM *sdm;
	Dfloat *BS0,*BN0,*BW0,*BE0,*BUp0,*BDown0;     // Boundaries Subdomains OUT
    Dfloat *BS1,*BN1,*BW1,*BE1,*BUp1,*BDown1;     // Boundaries Subdomains IN
    int N_SN,N_WE,N_UpDown;
    


public:
  
  MPI_DATA(SDM *Isdm);
  ~MPI_DATA();
  void TRANSFER(char *VarName);
 
};



#endif

#ifndef __FOURIER__
#define __FOURIER__

#include "definitions.hpp"
#include "sdm.hpp"

class DFT {

protected:
  int nf;
  SDM *sdm;

public:

  Dfloat **Fux,**Fuy,**Fuz;
  Dfloat **iFux,**iFuy,**iFuz;
  Dfloat *freq;

  DFT(SDM *in_sdm,std::string nFile,int infreq);
  ~DFT();

  void FD(Dfloat dt,int ktime);

  void InitVar(Dfloat f);
  
  
};


#endif

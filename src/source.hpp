// -*- c++ -*-


#ifndef __SOURCE__
#define __SOURCE__

#include <stdio.h>
#include "definitions.hpp"
#include "geometry3D.hpp"
#include <math.h>

class source {

protected:

  geometry3D *GDomain;
  std::string FileS,line;
  std::ifstream R;
  char data[200],*token;

public:

  Dfloat *Mxx,*Myy,*Mzz,*Mxy,*Mxz,*Myz;
  Dfloat *strike,*dip,*slip,*azimuth,*M0,*d_t0;
  Dfloat *xcoord,*ycoord,*zcoord;
  VecI *pos_src;
  int ns;

  source(geometry3D* domain, std::string nFile,int nsource);

  ~source();

  Dfloat sourceType(Dfloat t0,Dfloat f0, int itime,Dfloat dt, int Tsrc);

  
  void smoment();

  
  

};


#endif

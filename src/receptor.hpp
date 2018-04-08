// -*- c++ -*-


#ifndef __RECEPTOR__
#define __RECEPTOR__

#include <stdio.h>
#include "definitions.hpp"
#include "geometry3D.hpp"
#include <math.h>

class receptor {

protected:

  geometry3D *GDomain;
  std::string FileS,line;
  std::ifstream R;
  std::fstream *RX,*RY,*RZ;

public:

  Dfloat *xcoord,*ycoord,*zcoord;
  std::string *nameStation;
  VecI *pos_recep;
  int nr;
  Dfloat *vx_ad,*vy_ad,*vz_ad;
  int nt;
  

  receptor(geometry3D* domain, std::string nFile,int nrecep,int in_nt);

  void FileOpen(int i,int PROPAGATION);
  
  void FileClose(int i);

  void WriteFile(int i, Dfloat vx, Dfloat vy,Dfloat vz);

  void LoadFile(int i);

  ~receptor();
  
};


#endif

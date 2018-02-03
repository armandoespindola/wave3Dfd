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
  std::ofstream *RX,*RY,*RZ;

public:

  Dfloat *xcoord,*ycoord,*zcoord;
  std::string *nameStation;
  VecI *pos_recep;
  int nr;

  receptor(geometry3D* domain, std::string nFile,int nrecep);

  void FileOpen(int i);
  
  void FileClose(int i);

  void WriteFile(int i, Dfloat vx, Dfloat vy,Dfloat vz);

  ~receptor();
  
};


#endif

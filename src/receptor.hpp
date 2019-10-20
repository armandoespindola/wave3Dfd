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
  VecI *pos_vx;
  VecI *pos_vy;
  VecI *pos_vz;
  int nr;
  Dfloat *vx_ad,*vy_ad,*vz_ad;
  int nt;
  

  receptor(geometry3D* domain, std::string nFile,int nrecep,int in_nt);

  void FileOpen(int i,int PROPAGATION);
  
  void FileClose(int i);

  void WriteFile(int i, Dfloat *vx, Dfloat *vy,Dfloat *vz);

  void LoadFile(int i);

  ~receptor();

  // PRINT INFORMATION ABOUT RECEPTORS
  void PrintInf();
  
};


#endif

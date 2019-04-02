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


#ifndef __SOURCE__
#define __SOURCE__

#include <stdio.h>
#include "definitions.hpp"
#include "geometry3D.hpp"
#include "utilities.hpp"
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
  Dfloat sinc_wx[8],sinc_wy[8],sinc_wz[8];
  VecI *pos_src;
  int ns;
  int *nshift;

  source(geometry3D* domain, std::string nFile,int nsource);

  ~source();

  Dfloat sourceType(Dfloat t0,Dfloat f0, int itime,Dfloat dt, int Tsrc);

  
  void smoment();

  // PRINT INFORMATION ABOUT SOURCES
  void PrintInf();

  // WEIGHT SINC FUNCTION
  void w_sinc(int freeSurf,int is);

  
  

};


#endif

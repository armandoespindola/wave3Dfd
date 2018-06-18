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

#ifndef __FOURIER__
#define __FOURIER__

#include "definitions.hpp"
#include "sdm.hpp"

class DFT {

protected:
  int nf;
  SDM *sdm;

public:

  Dfloat **Fvx,**Fvy,**Fvz;
  Dfloat **iFvx,**iFvy,**iFvz;
  Dfloat *freq;

  DFT(SDM *in_sdm,std::string nFile,int infreq);
  ~DFT();

  void FD(Dfloat dt,int ktime);

  void InitVar(Dfloat f);
  
  
};


#endif

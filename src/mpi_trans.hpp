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

#ifndef __MPITRANS__
#define __MPITRANS__

#include "definitions.hpp"
#include "geometry3D.hpp"
#include "model.hpp"
#include "show.hpp"
#include "sdm.hpp"
#include "kernels.hpp"



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
  void Merge(Dfloat *LocVar,int NX,int NY, int NZ,VecI *subi,int rank);
  void MergePrint(Dfloat *LocVar,int NX,int NY, int NZ,VecI *subi,int rank, char *name);
};



#endif

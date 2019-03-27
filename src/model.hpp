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

  // CFL Condition Test
  int CFL(Dfloat dt,Dfloat dx,Dfloat dy,Dfloat dz);
  

};


#endif

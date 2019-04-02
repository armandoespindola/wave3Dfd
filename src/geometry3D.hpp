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


#ifndef __geometry3D__
#define __geometry3D__

#include "definitions.hpp"

// CLASS GEOMETRY 3D

class geometry3D {

protected:
  // GEOMETRY DIMENSION
  int dim;
  // NUMBER OF NODES
  VecI nod;
  // NUMBER OF ELEMENTS
  VecI ne;
  // TOTAL NUMBER OF ELEMENTS AND NODES
  VecI NE,NOD;
  // DOMAIN LIMITS
  VecF limI,limF;
  // DELTA SPACE x, y, z
  VecF Delta;
  // PML NODES
  VecI HALO;
  

public:

// CONSTRUCTOR WITHOUT PML
  geometry3D(VecF IlimI,VecF IlimF,VecI InElem);

  geometry3D(VecF IlimI,VecF IlimF,VecI InElem,VecI HALO);


  
// RETURN DIMENSION
  inline int Dim() { return dim;}

  // RETURN NUMBER OF TOTAL NODES
  inline int L_Node() { return (nod.x * nod.y * nod.z); }

  // NUMBER OF NODES IN X DIRECTION
  inline int L_NodeX() { return nod.x; }

  // NUMBER OF NODES IN Y DIRECTION
  inline int L_NodeY() { return nod.y; }

  // NUMBER OF NODES IN Z DIRECTION
  inline int L_NodeZ() { return nod.z; }

  // NUMBER OF ELEMENTS IN X DIRECTION
  inline int L_NelemX() { return ne.x; }

  // NUMBER OF ELEMENTS IN Y DIRECTION
  inline int L_NelemY() { return ne.y; }

  // NUMBER OF ELEMENTS IN Z DIRECTION

  inline int L_NelemZ() { return ne.z; } 

  // RETURN DELTA X
  inline Dfloat Dx() { return Delta.x;}

  // RETURN DELTA Y
  inline Dfloat Dy() { return Delta.y;}

  // RETURN DELTA Z
  inline Dfloat Dz() { return Delta.z;}


  inline int HALO_Nelem() { return (NE.x * NE.y * NE.z); }
  // NUMBER OF ELEMENTS WITH (HALO OR PML)
  inline int HALO_NelemX() { return NE.x;}

  // NUMBER OF ELEMENTS WITH (HALO OR PML)
  inline int HALO_NelemY() { return NE.y;}

  // NUMBER OF ELEMENTS WITH (HALO OR PML)
  inline int HALO_NelemZ() { return NE.z;}

  // RETURN NUMBER OF TOTAL NODES WITH (HALO OR PML)
  inline int HALO_Node() { return (NOD.x * NOD.y * NOD.z); }

  // NUMBER OF NODES IN X DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeX() { return NOD.x; }

  // NUMBER OF NODES IN Y DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeY() { return NOD.y; }

  // NUMBER OF NODES IN Z DIRECTION WITH (HALO OR PML)
  inline int HALO_NodeZ() { return NOD.z; }

  // X COORDINATE

  Dfloat CoorX(int indx);

  // Y COORDINATE

  Dfloat CoorY(int indx);

  // Z Coordinate

  Dfloat CoorZ(int indx);


  // X COORDINATE

  Dfloat CoorXHalf(int indx);

  // Y COORDINATE

  Dfloat CoorYHalf(int indx);

  // Z COORDINATE

  Dfloat CoorZHalf(int indx);

  // PML THICKNESS
  VecF thickness_PML();

  // HALO THICKNESS
  VecF thickness_HALO();
  

  // FIND NODE
  VecI FindNode(VecF coord);

  // FIND NODE GENERAL LOCATION
  VecI FindNode(VecF coord, VecI HALF);
  

  
  

  


  

};



#endif

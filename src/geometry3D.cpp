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


#include "geometry3D.hpp"
#include "definitions.hpp"



geometry3D::geometry3D(VecF IlimI,VecF IlimF,VecI InElem){
  
  // PML OR HALO NODES
  HALO.x = PML.x;
  HALO.y = PML.y;
  HALO.z = PML.z;

  // NUMBER ELEMENTS
  ne.x = InElem.x ;
  ne.y = InElem.y ;
  ne.z = InElem.z ;

  // NUMBER OF NODES ORIGINAL DOMAIN
  nod.x = InElem.x + 1;
  nod.y = InElem.y + 1;
  nod.z = InElem.z + 1;  

  // INITIAL LIMIT DOMAIN 
  limI.x = IlimI.x;
  limI.y = IlimI.y;
  limI.z = IlimI.z;

  // END LIMIT DOMAIN
  limF.x = IlimF.x;
  limF.y = IlimF.y;
  limF.z = IlimF.z;

  // DELTA DISCRETIZATION
  Delta.x = (limF.x - limI.x) / Dfloat(ne.x);
  Delta.y = (limF.y - limI.y) / Dfloat(ne.y);
  Delta.z = (limF.z - limI.z) / Dfloat(ne.z);

  // NUMBER ELEMENTS WITH PML OR HALO
  NE.x = InElem.x + 2 * HALO.x;
  NE.y = InElem.y + 2 * HALO.y;
  NE.z = InElem.z + HALO.z;

  // NUMBER OF NODES WITH PML OR HALO
  NOD.x = NE.x + 1;
  NOD.y = NE.y + 1;
  NOD.z = NE.z + 1;  
  
}



geometry3D::geometry3D(VecF IlimI,VecF IlimF,VecI InElem, VecI IHALO){
  

  // PML OR HALO NODES
  HALO.x = IHALO.x;
  HALO.y = IHALO.y;
  HALO.z = IHALO.z;

  // NUMBER ELEMENTS
  ne.x = InElem.x;
  ne.y = InElem.y;
  ne.z = InElem.z;

  // NUMBER OF NODES ORIGINAL DOMAIN
  nod.x = InElem.x + 1;
  nod.y = InElem.y + 1;
  nod.z = InElem.z + 1;  

  // INITIAL LIMIT DOMAIN 
  limI.x = IlimI.x;
  limI.y = IlimI.y;
  limI.z = IlimI.z;

  // END LIMIT DOMAIN
  limF.x = IlimF.x;
  limF.y = IlimF.y;
  limF.z = IlimF.z;

  // DELTA DISCRETIZATION
  Delta.x = (limF.x - limI.x) / Dfloat(ne.x);
  Delta.y = (limF.y - limI.y) / Dfloat(ne.y);
  Delta.z = (limF.z - limI.z) / Dfloat(ne.z);

// NUMBER ELEMENTS WITH PML OR HALO
  NE.x = InElem.x + 2 * HALO.x;
  NE.y = InElem.y + 2 * HALO.y;
  NE.z = InElem.z + 2 * HALO.z;

  // NUMBER OF NODES WITH PML OR HALO
  NOD.x = NE.x + 1;
  NOD.y = NE.y + 1;
  NOD.z = NE.z + 1;  

}






// COORDINATE ORIGINAL DOMAIN

Dfloat geometry3D::CoorX(int i){
  Dfloat coord;

  if (i < NOD.x){

    coord = limI.x + (Dfloat) i * Dx() - thickness_HALO().x;

  return coord;

  } else {
    
    std::cout<<"Out of Dimension_X"<<std::endl;

  }
}


Dfloat geometry3D::CoorY(int i){
  Dfloat coord;

  if (i < NOD.y){
    
    coord = limI.y + (Dfloat) i * Dy() - thickness_HALO().y;

  return coord;

  }else{
    
    std::cout<<"Out of Dimension_Y"<<std::endl;
    
  }   
}


Dfloat geometry3D::CoorZ(int i){
  Dfloat coord;

  if (i < NOD.z) {
    coord = limI.z + (Dfloat) i * Dz();
    return coord;
    
  }else{

    std::cout<<"Out of Dimension_Z"<<std::endl;

  }
}



Dfloat geometry3D::CoorXHalf(int i){
  Dfloat coord;

  if (i < NOD.x-1){

    coord = limI.x + (Dx() / 2.0) + (Dfloat) i * Dx() - thickness_HALO().x;

  return coord;

  } else {
    
    std::cout<<"Out of Dimension_X_half"<<std::endl;

  }  
}


Dfloat geometry3D::CoorYHalf(int i){
  Dfloat coord;

  if (i < NOD.y-1){
    
    coord = limI.y + (Dy() / 2.0) + (Dfloat) i * Dy() - thickness_HALO().y;

  return coord;

  }else{
    
    std::cout<<"Out of Dimension_Y_half"<<std::endl;
    
  }    
}


Dfloat geometry3D::CoorZHalf(int i){
  Dfloat coord;

  if (i < NOD.z-1) {
    coord = limI.z - (Dz() / 2.0) + (Dfloat) i * Dz();
    return coord;
    
  }else{

    std::cout<<"Out of Dimension_Z_half"<<std::endl;
  }
}


// THICKNES PML 

VecF geometry3D::thickness_PML(){

  VecF thickness;

  thickness.x =  (Dfloat) PML.x * Delta.x;
  thickness.y =  (Dfloat) PML.y * Delta.y;
  thickness.z =  (Dfloat) PML.z * Delta.z;

  return thickness;
  
}

// THICKNES HALO

VecF geometry3D::thickness_HALO(){

  VecF thickness;

  thickness.x =  (Dfloat) HALO.x * Delta.x;
  thickness.y =  (Dfloat) HALO.y * Delta.y;
  thickness.z =  (Dfloat) HALO.z * Delta.z;

  return thickness;
  
}


VecI geometry3D::FindNode(VecF coord){

  Dfloat deltaAvg = (Delta.x +  Delta.z + Delta.y) / 3.0;
  Dfloat dist2,dist1;
  VecI ind;


dist1 = deltaAvg;
for (int iz=0;iz<NOD.z;iz++){
if ( (CoorZ(iz) >= limI.z) && (CoorZ(iz) <= limF.z)){
dist2 = abs(coord.z - CoorZ(iz));
if (dist2<dist1){
dist1 = dist2;
ind.z = (NOD.z - 1) - iz;
}
}
}



dist1 = deltaAvg;
for (int ix=0;ix<NOD.x;ix++){
if ( (CoorX(ix) >= limI.x) && (CoorX(ix) <= limF.x)){
dist2 = abs(coord.x - CoorX(ix));
if (dist2<dist1){
dist1 = dist2;
ind.x = ix;
}
}
}

dist1 = deltaAvg;
for (int iy=0;iy<NOD.y;iy++){
if ( (CoorY(iy) >= limI.y) && (CoorY(iy) <= limF.y)){
dist2 = abs(coord.y - CoorY(iy));
if (dist2<dist1){
dist1 = dist2;
ind.y = iy;
}
}
}

return ind;

}


VecI geometry3D::FindNode(VecF coord, VecI HALF){

  Dfloat deltaAvg = (Delta.x +  Delta.z + Delta.y) / 3.0;
  Dfloat dist2,dist1;
  VecI ind;


  if (HALF.z){
dist1 = deltaAvg;
for (int iz=0;iz<NOD.z-1;iz++){
if ( (CoorZHalf(iz) >= limI.z) && (CoorZHalf(iz) <= limF.z)){
dist2 = abs(coord.z - CoorZHalf(iz));
if (dist2<dist1){
dist1 = dist2;
ind.z = (NOD.z - 1) - iz;
}
}
}

  } else {

dist1 = deltaAvg;
for (int iz=0;iz<NOD.z;iz++){
if ( (CoorZ(iz) >= limI.z) && (CoorZ(iz) <= limF.z)){
dist2 = abs(coord.z - CoorZ(iz));
if (dist2<dist1){
dist1 = dist2;
ind.z = (NOD.z - 1) - iz;
}
}
}

  }

  
  if (HALF.x){
dist1 = deltaAvg;
for (int ix=0;ix<NOD.x-1;ix++){
if ( (CoorXHalf(ix) >= limI.x) && (CoorXHalf(ix) <= limF.x)){
dist2 = abs(coord.x - CoorXHalf(ix));
if (dist2<dist1){
dist1 = dist2;
ind.x = ix;
}
}
}

  } else {


dist1 = deltaAvg;
for (int ix=0;ix<NOD.x;ix++){
if ( (CoorX(ix) >= limI.x) && (CoorX(ix) <= limF.x)){
dist2 = abs(coord.x - CoorX(ix));
if (dist2<dist1){
dist1 = dist2;
ind.x = ix;
}
}
}
 

  }


  if (HALF.y) {
dist1 = deltaAvg;
for (int iy=0;iy<NOD.y-1;iy++){
if ( (CoorYHalf(iy) >= limI.y) && (CoorYHalf(iy) <= limF.y)){
dist2 = abs(coord.y - CoorYHalf(iy));
if (dist2<dist1){
dist1 = dist2;
ind.y = iy;
}
}
}
  } else {



dist1 = deltaAvg;
for (int iy=0;iy<NOD.y;iy++){
if ( (CoorY(iy) >= limI.y) && (CoorY(iy) <= limF.y)){
dist2 = abs(coord.y - CoorY(iy));
if (dist2<dist1){
dist1 = dist2;
ind.y = iy;
}
}
}
 

  }

  
return ind;

}




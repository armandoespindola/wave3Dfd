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

 #include "definitions.hpp"
 #include "pml.hpp"

  // NX NODES WITHOUT PML
  // PML_NX NODES WITH PML
  // T_PML  THICKNESS PML
  // Delta  DELTA
  // DPML   NODES PML
  // dt   DELTA TIME
  // f0   SOURCE FREQUENCY
  // bool MIN_PML
  // bool MAX_PML

DPML::DPML(int NX,int PML_NX, Dfloat T_PML, Dfloat Delta, int IPML,Dfloat dt, \
 Dfloat f0, bool PML_XMIN_USE, bool PML_XMAX_USE){

DELTAX = Delta;
NXT = PML_NX;
NPOWER = 2.0;
cp = 7000.0;
Rcoef = 0.001;
ALPHA_MAX_PML = 2.0 * pi * (f0 / 2.0);
n_pml = IPML;
n_x = NX;
thickness_PML_x = T_PML;
// # DAMPING PROFILE  
K_x  = ONE;
d0_x= - 1.0 * (NPOWER +1.0) * cp * log(Rcoef) / (2.0 * thickness_PML_x);

  // # ORIGIN OF THE PML LAYER
X_OrgL = thickness_PML_x;
X_OrgR = (Dfloat)(NXT - 1.0) * DELTAX - thickness_PML_x;
X_OrgR = (Dfloat)(n_x+n_pml-1)*DELTAX;


d_x = new Dfloat [PML_NX]; 
d_x_half = new Dfloat [PML_NX];
alpha_x = new Dfloat [PML_NX];
alpha_x_half = new Dfloat [PML_NX];
a_x = new Dfloat [PML_NX];
a_x_half = new Dfloat [PML_NX];
b_x = new Dfloat [PML_NX];
b_x_half = new Dfloat [PML_NX];



for(int i=0;i<NXT;i++){
  d_x[i] = ZERO; 
  d_x_half[i] = ZERO;
  alpha_x[i] = ZERO;
  alpha_x_half[i] = ZERO;
  a_x[i] = ZERO;
  a_x_half[i] = ZERO;
  b_x[i] = ZERO;
  b_x_half[i] = ZERO;
}


for (int i=0;i<NXT;i++){

  
  // # LEFT EDGE PML X
    if ( PML_XMIN_USE) {

  // Abscissa of current grid point along the damping profile
  val = X_OrgL - (DELTAX * (Dfloat)i);

    // DAMPING PROFLE AT THE GRID POINTS
    if (val >= ZERO)  {
  val_NORM = val / thickness_PML_x;
  d_x[i] = d0_x * pow(val_NORM,NPOWER);
  alpha_x[i] = ALPHA_MAX_PML * (ONE - val_NORM) + 0.10 * ALPHA_MAX_PML;
  b_x[i] = exp(-1.0 * (d_x[i] / K_x + alpha_x[i]) * dt);
  a_x[i] = (d_x[i] / ( K_x * (d_x[i] + K_x * alpha_x[i]) ) ) * (b_x[i] - 1.0);

  //printf("%i\t%f\t%f\t%f\t%f\n",i,d_x[i],alpha_x[i],b_x[i],a_x[i]);

  
    }

    
  // # DAMPING PROFLE AT HALF THE GRID POINTS

  // Abscissa of current grid point along the damping profile
  val = X_OrgL - (DELTAX * (Dfloat)i + DELTAX / 2.0 );

  // #  DAMPING PROFLE AT THE GRID POINTS
    if (val >= ZERO)  {
  val_NORM = val / thickness_PML_x;
  d_x_half[i] = d0_x * pow(val_NORM,NPOWER);
  alpha_x_half[i] = ALPHA_MAX_PML * (ONE - val_NORM) + 0.10 * ALPHA_MAX_PML;
  b_x_half[i] = exp(-1.0 * (d_x_half[i] / K_x + alpha_x_half[i]) * dt);
  a_x_half[i] = (d_x_half[i] / ( K_x * (d_x_half[i] + K_x * alpha_x_half[i]) ) ) * (b_x_half[i] - 1.0);
  
  //printf("Mid-Points %i\t%f\t%f\t%f\t%f\n",i,d_x_half[i],alpha_x_half[i],b_x_half[i],a_x_half[i]);
}
  
    }
    
 if ( PML_XMAX_USE) {


   // Abscissa of current grid point along the damping profile
  val = (DELTAX * (Dfloat)i) - X_OrgR;

    // DAMPING PROFLE AT THE GRID POINTS
    if (val >= ZERO)  {
  val_NORM = val / thickness_PML_x;
  d_x[i] = d0_x * pow(val_NORM,NPOWER);
  alpha_x[i] = ALPHA_MAX_PML * (ONE - val_NORM) + 0.10 * ALPHA_MAX_PML;
  b_x[i] = exp(-1.0 * (d_x[i] / K_x + alpha_x[i]) * dt);
  a_x[i] = (d_x[i] / ( K_x * (d_x[i] + K_x * alpha_x[i]) ) ) * (b_x[i] - 1.0);

  //  printf("%i\t%f\t%f\t%f\t%f\n",i,d_x[i],alpha_x[i],b_x[i],a_x[i]); 
}

  // # DAMPING PROFLE AT HALF THE GRID POINTS

  // Abscissa of current grid point along the damping profile
  val = (DELTAX * (Dfloat)i + DELTAX / 2.0 ) - X_OrgR;

  // #  DAMPING PROFLE AT THE GRID POINTS
    if (val > 0.0)  {
  val_NORM = val / thickness_PML_x;
  d_x_half[i] = d0_x * pow(val_NORM,NPOWER);
  alpha_x_half[i] = ALPHA_MAX_PML * (ONE - val_NORM) + 0.10 * ALPHA_MAX_PML;
  b_x_half[i] = exp(-1.0 * (d_x_half[i] / K_x + alpha_x_half[i]) * dt);
  a_x_half[i] = (d_x_half[i] / ( K_x * (d_x_half[i] + K_x * alpha_x_half[i]) ) ) * (b_x_half[i] - 1.0);
  
  // printf("Mid-Points %i\t%f\t%f\t%f\t%f\n",i,d_x_half[i],alpha_x_half[i],b_x_half[i],a_x_half[i]);
}
 }

 // printf("%i\t%f\t%f\t%f\t%f\n",i,d_x[i],alpha_x[i],b_x[i],a_x[i]);
 // printf("Mid-Points %i\t%f\t%f\t%f\t%f\n",i,d_x_half[i],alpha_x_half[i],b_x_half[i],a_x_half[i]);
    }


}


DPML::~DPML(){


delete [] d_x;
delete [] d_x_half;
delete [] alpha_x;
delete [] alpha_x_half;
delete [] a_x;
delete [] a_x_half;
delete [] b_x;
delete [] b_x_half;

}


Dfloat DPML::B(int i){

  if (i < NXT) {

    return b_x[i];

  }else{

    std::cout<<"Out of Dimension_PML"<<std::endl;
  }


}

Dfloat DPML::B_HALF(int i){

  if (i < NXT) {

    return b_x_half[i];

  }else{

    std::cout<<"Out of Dimension_PML"<<std::endl;
  }


}

Dfloat DPML::A(int i){

  if (i < NXT) {

    return a_x[i];

  }else{

    std::cout<<"Out of Dimension_PML"<<std::endl;
  }


}

Dfloat DPML::A_HALF(int i){

  if (i < NXT) {

    return a_x_half[i];

  }else{

    std::cout<<"Out of Dimension_PML"<<std::endl;
  }


}






 

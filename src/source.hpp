// -*- c++ -*-
#ifndef __SOURCE__
#define __SOURCE__
#include "source.hpp"

#include <stdio.h>
#include "definitions.hpp"
#include "parameters.hpp"
#include <math.h>


Dfloat source(Dfloat t0 , int itime) {

  Dfloat src,a_fu,amp,time;

  time = 0.5 * dt + itime  * dt; 
  a_fu= pow (pi*f0,2.0);
  src = 0.0;

  // GAUSSIAN 
  
  if (T_SRC==0){ 
    src = exp(-a_fu * pow(time - t0,2.0));
  }
  
  // FIRST DERIVATIVE OF A GAUSSIAN
  
  if (T_SRC==1){
    src = 4.0 * a_fu *(time - t0) * exp(-2.0 * a_fu * pow( (time  - t0),2.0) );
  }

  // SECOND DERIVATIVE OF A GAUSSIAN (RICKER PULSE)
  
  if (T_SRC==2){
    src = (1.0 - 2.0 * a_fu * pow((time - t0),2.0)) * exp(-a_fu * pow((time - t0),2.0)) ;
  }

  // HEAVISIDE STEP FUNCTION
  
  if (T_SRC==3){
    src = 1.0 * 1.0e-5 * (1.0/100.0); //dyn*cm -> N*m
  }

  
  return src;

};

#endif

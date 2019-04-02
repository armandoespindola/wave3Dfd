/* bessel.c
                      Copyright (c) 1998 
                  Kapteyn Institute Groningen
                     All Rights Reserved.
*/

/*
#>            bessel.dc2

Function:     BESSEL

Purpose:      Evaluate Bessel function J, Y, I, K of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          See bessj.dc2, bessy.dc2, bessi.dc2 or bessk.dc2
                      
Description:  The differential equation 

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- + (x - v ).w = 0
                        2       dx
                      dx
                      
              has two solutions called Bessel functions of the first kind
              Jv(x) and Bessel functions of the second kind Yv(x).
              The routines bessj and bessy return the J and Y for 
              integer v and therefore are called Bessel functions 
              of integer order.
              
              The differential equation 

                       2
                   2  d w       dw      2   2
                  x . --- + x . --- - (x + v ).w = 0
                        2       dx
                      dx
                      
              has two solutions called modified Bessel functions
              Iv(x) and Kv(x).
              The routines bessi and bessk return the I and K for 
              integer v and therefore are called Modified Bessel 
              functions of integer order.
              (Abramowitz & Stegun, Handbook of mathematical 
              functions, ch. 9, pages 358,- and 374,- )
                            
              The implementation is based on the ideas from 
              Numerical Recipes, Press et. al. 
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
*/

#ifndef __UTILITIES__
#define __UTILITIES__

#include <stdio.h>
#include "definitions.hpp"
#include <math.h>

Dfloat bessj0(Dfloat x);
Dfloat bessi0( Dfloat x );

#endif

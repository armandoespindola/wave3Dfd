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

#ifndef __definitions__
#define __definitions__


// PRECISION OF OPERATIONS //

#define ZERO 0.0
#define ONE  1.0
#define DREAL float
typedef DREAL Dfloat;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <string>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <iomanip>
#include <fstream>

#define MY_MPI_Dfloat MPI_FLOAT

//  VEC 3 ELEMENTS FLOAT
struct VecF{
	Dfloat x;
	Dfloat y;
	Dfloat z;
};


//  VEC 3 ELEMENTS INT
struct VecI{
	int x;
	int y;
	int z;
};


// PML PARAMETERS

const bool PML_XMIN = true;
const bool PML_XMAX = true;
const bool PML_YMIN = true;
const bool PML_YMAX = true;
const bool PML_ZMIN = true;
const bool PML_ZMAX = false;

// PML NODES

//const VecI PML={16,16,16};


// # FINITE DIFFERENCES PARAMETERS 4TH ORDER SPACE

const int KHALO = 2;  // Half of the order on space (k/2)

// STENCIL COEFFICIENTS

const Dfloat C0 = 1.0 / 24.0;
const Dfloat C1 = 9.0 / 8.0;

// CONSTANTS

const Dfloat pi = 3.141592653589793;


// SAVE BOUNDARIES INTERVAL

const int stepb = 100;


// MACROS LOOPS

#define FOR_IJK(iz,zinit,zend,iy,yinit,yend,ix,xinit,xend) for(int iz=((int)zinit);iz<((int)zend);++iz){ \
  for(int iy=((int)yinit);iy<((int)yend);++iy){ for(int ix=((int)xinit);ix<((int)xend);++ix){;

#define END_FOR_IJK } } }

#define FOR_IJ(iy,yinit,yend,ix,xinit,xend) for(int iy=((int)yinit);iy<((int)yend);++iy){ \
  for(int ix=((int)xinit);ix<((int)xend);++ix){

#define END_FOR_IJ } }
#endif

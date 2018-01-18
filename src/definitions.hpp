#ifndef __definitions__
#define __definitions__


// PRECISION OF OPERATIONS //

#define ZERO 0.0
#define ONE  1.0
#define DREAL double
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

const VecF PML={16,16,16};


// # FINITE DIFFERENCES PARAMETERS 4TH ORDER SPACE

const int KHALO = 2;  // Half of the order on space (k/2)

// STENCIL COEFFICIENTS

const Dfloat C0 = 1.0 / 24.0;
const Dfloat C1 = 9.0 / 8.0;

// CONSTANTS

const Dfloat pi = 3.141592653589793;



#endif

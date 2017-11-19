#ifndef __SDM__
#define __SDM__

#include "definitions.hpp"
#include "pml.hpp"
#include "geometry3D.hpp"

class SDM {

protected:


        // GI INITIAL GLOBAL LIMIT
	// GF END GLOBAL LIMIT 
	// NodG NUMBER OF GLOBAL NODES
	// IlimI INITIAL LIMIT LOCAL DOMAIN
	// IlimF END LIMIT LOCAL DOMAIN
	// INOD   NUMBER OF LOCAL NODES
	// f0 FREQUENCY
	// dt DELTA T
        // Nsdm SUBDOMAIN INDEX

	VecF GI,GF;
	DPML *pml_x,*pml_y,*pml_z;
	geometry3D *SDMGeom;
	VecI HALO,NodG;
	Dfloat f0,dt;
	VecF thickness_PML;
	VecI NT,Nsdm;
	VecI NodLoc;

	// WAVE PROPAGATION VARIABLES  
	Dfloat *sxx,*syy,*szz;
	Dfloat *sxy,*sxz,*syz;
	Dfloat *vx,*vy,*vz;
	Dfloat *ux,*uy,*uz;
	Dfloat *mu,*lamb,*rho;

	 // ### SHARED VARIBALES PML DAMPING
	Dfloat *dsxx_dx,*dsxy_dy,*dsxz_dz;
	Dfloat *dsxy_dx,*dsyy_dy,*dsyz_dz;
	Dfloat *dsxz_dx,*dsyz_dy,*dszz_dz;

	Dfloat *dvx_dx,*dvy_dy,*dvz_dz;
	Dfloat *dvx_dy,*dvy_dx,*dvx_dz;
	Dfloat *dvz_dx,*dvy_dz,*dvz_dy;

public:

	// GI INITIAL GLOBAL LIMIT
	// GF END GLOBAL LIMIT 
	// NodG NUMBER OF GLOBAL NODES
	// IlimI INITIAL LIMIT LOCAL DOMAIN
	// IlimF END LIMIT LOCAL DOMAIN
	// INOD   NUMBER OF LOCAL NODES
	// f0 FREQUENCY
	// dt DELTA T 


  SDM(VecF IGI, VecF IGF,VecI I_NodG,VecF IlimI, VecF IlimF, VecI INod, \
		Dfloat If0, Dfloat Idt, VecI INsdm);


  ~SDM();

  // INDEX_IJK
  
  int IJK(int i,int j, int k);

  // MODEL READ MODEL

  void ModelRead(Dfloat *model, char *param);

  // STABILITY CONDITION
  int CFL();

  


  





};

#endif

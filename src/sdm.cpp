#include "sdm.hpp"

SDM::SDM(VecF IGI, VecF IGF,VecI I_NodG,VecF IlimI, VecF IlimF, VecI INod, \
		Dfloat If0, Dfloat Idt, VecI INsdm) {

	// GI INITIAL GLOBAL LIMIT
	// GF END GLOBAL LIMIT 
	// NodG NUMBER OF GLOBAL NODES
	// IlimI INITIAL LIMIT LOCAL DOMAIN
	// IlimF END LIMIT LOCAL DOMAIN
	// INOD   NUMBER OF LOCAL NODES
	// f0 FREQUENCY
	// dt DELTA T
    // Nsdm SUBDOMAIN INDEX

	GI = IGI;
	GF = IGF;
	f0 = If0;
	dt = Idt;
	NodG = I_NodG;
	HALO.x = KHALO;
	HALO.y = KHALO;
	HALO.z = KHALO;
	NodLoc = INod;  

	VecI ELE={INod.x-1,INod.y-1,INod.z-1};

	SDMGeom = new geometry3D(IlimI,IlimF,ELE,HALO);

	thickness_PML = SDMGeom->thickness_PML();

	NT.x = NodG.x + 2 * PML.x;
	NT.y = NodG.y + 2 * PML.y;
	NT.z = NodG.z + PML.z;

	pml_x = new DPML(NodG.x,NT.x,thickness_PML.x,SDMGeom->Dx(),PML.x,dt,f0, \
		PML_XMIN,PML_XMAX);

	pml_y = new DPML(NodG.y,NT.y,thickness_PML.y,SDMGeom->Dy(),PML.y,dt,f0, \
		PML_YMIN,PML_YMAX);

	pml_z = new DPML(NodG.z,NT.z,thickness_PML.z,SDMGeom->Dz(),PML.z,dt,f0, \
		PML_ZMIN,PML_ZMAX);



	sxx = new Dfloat [SDMGeom->HALO_Node()];
	syy = new Dfloat [SDMGeom->HALO_Node()];
	szz = new Dfloat [SDMGeom->HALO_Node()];
	sxy = new Dfloat [SDMGeom->HALO_Node()];
	sxz = new Dfloat [SDMGeom->HALO_Node()];
	syz = new Dfloat [SDMGeom->HALO_Node()];
	vx = new Dfloat [SDMGeom->HALO_Node()];
	vy = new Dfloat [SDMGeom->HALO_Node()];
	vz = new Dfloat [SDMGeom->HALO_Node()];
	ux = new Dfloat [SDMGeom->HALO_Node()];
	uy = new Dfloat [SDMGeom->HALO_Node()];
	uz = new Dfloat [SDMGeom->HALO_Node()];
	mu = new Dfloat [SDMGeom->HALO_Node()];
	lamb = new Dfloat [SDMGeom->HALO_Node()];
	rho = new Dfloat [SDMGeom->HALO_Node()];
	dsxx_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsxy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dsxz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dsxy_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsyy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dsyz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dsxz_dx = new Dfloat [SDMGeom->HALO_Node()];
	dsyz_dy = new Dfloat [SDMGeom->HALO_Node()];
	dszz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dy = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dy = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvx_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dx = new Dfloat [SDMGeom->HALO_Node()];
	dvy_dz = new Dfloat [SDMGeom->HALO_Node()];
	dvz_dy = new Dfloat [SDMGeom->HALO_Node()];


}


SDM::~SDM(){

	delete [] pml_x;
	delete [] pml_y;
	delete [] pml_z;
	delete [] SDMGeom;
	delete [] sxx;
	delete [] syy;
	delete [] szz;
	delete [] sxy;
	delete [] sxz;
	delete [] syz;
	delete [] vx;
	delete [] vy;
	delete [] vz;
	delete [] ux;
	delete [] uy;
	delete [] uz;
	delete [] mu;
	delete [] lamb;
	delete [] rho;
	delete [] dsxx_dx;
	delete [] dsxy_dy;
	delete [] dsxz_dz;
	delete [] dsxy_dx;
	delete [] dsyy_dy;
	delete [] dsyz_dz;
	delete [] dsxz_dx;
	delete [] dsyz_dy;
	delete [] dszz_dz;
	delete [] dvx_dx;
	delete [] dvy_dy;
	delete [] dvz_dz;
	delete [] dvx_dy;
	delete [] dvy_dx;
	delete [] dvx_dz;
	delete [] dvz_dx;
	delete [] dvy_dz;
	delete [] dvz_dy;

}


// INDEX_IJK

int SDM::IJK(int i, int j, int k){

  int indx = i + j * NodLoc.x + k * NodLoc.x * NodLoc.y;

  return indx;

}


void SDM::ModelRead(Dfloat *model,char *param){

  if (strcmp("RHO",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i)
      rho[i] = model[i];

  }

  
  if (strcmp("MU",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i)
      mu[i] = model[i];

  }

  
  if (strcmp("LAMB",param) == 0){

    for (int i=0; i<SDMGeom->HALO_Node(); ++i)
      lamb[i] = model[i];

  }

}


int SDM::CFL(){

	for (int i=0; i<SDMGeom->HALO_Node(); ++i){

		if ((mu[i] < 0.0) || (rho[i] < 0.0) || (lamb[i] < 0.0))
		{
			std::cout<<"ERROR IN MODEL PARAMETERS"<<std::endl;

		} else {


		       if ((mu[i] == 0.0) || (rho[i] == 0.0) || (lamb[i] == 0.0))
			   break;

			Dfloat vp;

			vp = sqrt((lamb[i] + 2.0 * mu[i]) / rho[i]);

			VecF K;

			K.x = SDMGeom->Dx() / ( dt * sqrt(2.0 * vp * (C0 + C1)));
			K.y = SDMGeom->Dy() / ( dt * sqrt(2.0 * vp * (C0 + C1)));
			K.z = SDMGeom->Dz() / ( dt * sqrt(2.0 * vp * (C0 + C1)));

			if ((K.x >= 1.0) || (K.y >= 1.0) || (K.z >= 1.0))
			{
				std::cout<<"CFL NOT SATISFIED"<<std::endl;
				return 0;	
			} else{

				return 1;
			}

		}


	}


}


#ifndef __SDM__
#define __SDM__

#include "definitions.hpp"
#include "pml.hpp"
#include "geometry3D.hpp"

class SDM {

protected:

        // GI INITIAL GLOBAL LIMIT WITH PML 
	// GF END GLOBAL LIMIT WITH PML
	// NodG NUMBER OF GLOBAL NODES WITHOUT PML
	// IlimI INITIAL LIMIT SUB DOMAIN 
	// IlimF END LIMIT SUB DOMAIN
	// INOD   NUMBER OF LOCAL NODES WITHOUT HALO SUBDOMAIN
	// f0 FREQUENCY
	// dt DELTA T
        // Nsdm SUBDOMAIN INDEX
        // TOTAL NUMBER SUBDOMAIN NumSubDom

	VecF GI,GF;
	DPML *pml_x,*pml_y,*pml_z;
	geometry3D *SDMGeom;
	VecI HALO,NodG;
	Dfloat f0,dt;
	VecF thickness_PML;
        VecI NT,Nsdm,NumSubDom;
	VecI NodLoc;
  int N_omp;

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
      Dfloat If0, Dfloat Idt, VecI INsdm,VecI INumSubDom);


  ~SDM();

  // INDEX_IJK
  
  int IJK(int i,int j, int k);

  // MODEL READ MODEL

  void ModelRead(Dfloat *model, char param[]);

  // STABILITY CONDITION
  int CFL();


  // Initialize Veriables Propagation

  void InitVar(Dfloat f);

  // EXPORT BOUNDARY
  void EB(Dfloat *BN, Dfloat *DomLoc, char *TBound);
  void ExpBoundary(Dfloat *BN, char *TBound,char *NameVar);

  // IMPORT BOUNDARY
  void IB(Dfloat *BN, Dfloat *DomLoc, char *TBound);
  void ImpBoundary(Dfloat *BN, char *TBound,char *NameVar);
  
   // X COORDINATE

  Dfloat SCoorX(int indx);

  // Y COORDINATE

  Dfloat SCoorY(int indx);

  // Z Coordinate

  Dfloat SCoorZ(int indx);

   // X COORDINATE

  Dfloat SCoorXHalf(int indx);

  // Y COORDINATE

  Dfloat SCoorYHalf(int indx);

  // Z COORDINATE

  Dfloat SCoorZHalf(int indx);

  // LOCAL NODE TO GLOBAL NODE

  VecI Loc2Glo(VecI indx);

  // FIND NODE IN DOMAIN  
  VecI SFindNode(VecI coord);

  // FIND NODE IN DOMAIN 
  //VecI SFindNodeHalf(VecI coord);

  // ADD VALUE 
  void AddVal(VecI Indx,char *NameVar, Dfloat Val);

  // GET VALUE
  Dfloat GetVal(VecI Indx,char *NameVar);

    // FINITE DIFFERENCE STRESS SII

  void FD_SII(VecI Init,VecI Iend);


  // FINITE DIFFERENCE VELOCITY VX,VY,VZ
  
  void FD_VX(VecI Init,VecI Iend);

  void FD_VY(VecI Init,VecI Iend);

  void FD_VZ(VecI Init,VecI Iend);

  void FD_SXY(VecI Init,VecI Iend);

  void FD_SXZ(VecI Init,VecI Iend);

  void FD_SYZ(VecI Init,VecI Iend);


  // FINITE DIFFERENCES SXX, SYY, SZZ
  void FDSII();

  void FreeS_SII(VecI Init ,VecI Iend);

  // FINITE DIFFERENCES SXY
  void FDSXY();

  // FINITE DIFFERENCES SXZ
  void FDSXZ();

  // FINITE DIFFERENCES SYZ
  void FDSYZ();

  // FINITE DIFFERENCES SVX
  void FDVX();

  // FINITE DIFFERENCES SVY
  void FDVY();

  // FINITE DIFFERENCES SVZ
  void FDVZ();

  // Source Function 
  Dfloat source(int T_SRC,int itime, Dfloat t0);

  void printfile(Dfloat * Var,char *nfile, int ktime);

  void print(char* NameVar,int time);

  int BNorth();
  int BSouth();
  int BEast();
  int BWest();
  int BDown();
  int BUp();

  int SNodeT(){

    return SDMGeom->HALO_Node();

  };

  // NUMBER OF NODES IN X DIRECTION WITH (HALO OR PML)
  inline int SNodeX() { return SDMGeom->HALO_NodeX(); }

  // NUMBER OF NODES IN Y DIRECTION WITH (HALO OR PML)
  inline int SNodeY() { return SDMGeom->HALO_NodeY(); }

  // NUMBER OF NODES IN Z DIRECTION WITH (HALO OR PML)
  inline int SNodeZ() { return SDMGeom->HALO_NodeZ(); }

  // NUMBER OF NODES IN X DIRECTION 
  inline int L_SNodeX() { return SDMGeom->L_NodeX(); }

  // NUMBER OF NODES IN Y DIRECTION 
  inline int L_SNodeY() { return SDMGeom->L_NodeY(); }

  // NUMBER OF NODES IN Z DIRECTION
  inline int L_SNodeZ() { return SDMGeom->L_NodeZ(); }

  // INDEX SUBDOMAIN
  inline VecI SubIdx() { return Nsdm; }

  // NUMBER SUBDOMAINS
  inline VecI SubNum() { return NumSubDom; }

  inline void set_omp(int num) {N_omp = num;}

  


  

  





  





};

#endif

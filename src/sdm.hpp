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

#ifndef __SDM__
#define __SDM__

#include "definitions.hpp"
#include "pml.hpp"
#include "geometry3D.hpp"
#include "source.hpp"
#include "receptor.hpp"

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
	// PROPAGATION; 0 FORWARD PROPAGATION; 1 REVERSE PROPAGATION

	VecF GI,GF;
	DPML *pml_x,*pml_y,*pml_z;
	VecI HALO,NodG;
	VecF thickness_PML;
        VecI NT;
	int PROPAGATION;
        Dfloat sgn;

	
	 // ### SHARED VARIBALES PML DAMPING
	Dfloat *dsxx_dx,*dsxy_dy,*dsxz_dz;
	Dfloat *dsxy_dx,*dsyy_dy,*dsyz_dz;
	Dfloat *dsxz_dx,*dsyz_dy,*dszz_dz;

	Dfloat *dvx_dx,*dvy_dy,*dvz_dz;
	Dfloat *dvx_dy,*dvy_dx,*dvx_dz;
	Dfloat *dvz_dx,*dvy_dz,*dvz_dy;
        Dfloat *Rvx,*Rvy,*Rvz;

  // SOURCE MOMENT TENSOR

  source *sourceM;
  VecI *idx_source;
  int ntimesrc;
  // RECEPTORS

  receptor *station;
  VecI *idx_station;

  // BOUNDARIES RETRO-PROPAGATION
  Dfloat *bn_lx,*bn_rx;
  Dfloat *bn_ly,*bn_ry;
  Dfloat *bn_lz;
  // SAVE BOUNDARIES N-TIME STEPS
  int nsteps = stepb;

  // MPI FILES

  MPI_Status status;
  MPI_File BLX,BRX,BLY,BRY,BLZ;

  

public:

  // GEOMETRY MODEL
  geometry3D *SDMGeom;
  Dfloat f0,dt;
  int N_omp;
  
  // WAVE PROPAGATION VARIABLES  
	Dfloat *sxx,*syy,*szz;
	Dfloat *sxy,*sxz,*syz;
	Dfloat *vx,*vy,*vz;
	Dfloat *ux,*uy,*uz;
	Dfloat *mu,*lamb,*rho;

  VecI NumSubDom,Nsdm;
	VecI NodLoc;
  std::string FileSrc;
  int FileSrcB;
  Dfloat *srct;

	// GI INITIAL GLOBAL LIMIT
	// GF END GLOBAL LIMIT 
	// NodG NUMBER OF GLOBAL NODES
	// IlimI INITIAL LIMIT LOCAL DOMAIN
	// IlimF END LIMIT LOCAL DOMAIN
	// INOD   NUMBER OF LOCAL NODES
	// f0 FREQUENCY
	// dt DELTA T 


  SDM(VecF IGI, VecF IGF,VecI I_NodG,VecF IlimI, VecF IlimF, VecI INod, \
      Dfloat If0, Dfloat Idt, VecI INsdm,VecI INumSubDom,int iPROPAGATION);


  ~SDM();

  // INDEX_IJK
  
  int IJK(int i,int j, int k);

  // MODEL READ MODEL

  void ModelRead(Dfloat *model, char param[]);

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
  void Free_SII(VecI Init ,VecI Iend, int );
  void FDSII();

  // FINITE DIFFERENCES SXY
  void FDSXY();

  // FINITE DIFFERENCES SXZ
  void Free_SXZ(VecI Init ,VecI Iend);
  void FDSXZ();

  // FINITE DIFFERENCES SYZ
  void Free_SYZ(VecI Init ,VecI Iend);
  void FDSYZ();

  // FINITE DIFFERENCES SVX
  void Free_VX(VecI Init ,VecI Iend, int);
  void FDVX();

  // FINITE DIFFERENCES SVY
  void Free_VY(VecI Init ,VecI Iend, int);
  void FDVY();

  // FINITE DIFFERENCES SVZ
  void Free_VZ(VecI Init ,VecI Iend);
  void FDVZ();

  // SOURCE MOMENT TENSOR

  // INITIALIZE SOURCE
  void InitSource(geometry3D *GDomain,std::string nFile,int nsource,int SrcFile,int nt);

  // ADD VALUES SOURCE
  void AddSource(int itime,int T_SRC);

  // FINALIZE SOURCES

  void EndSource();

  // ADD AJOINT SOURCE

  void AddSourceAdj(int itime);

  // INITIALIZE RECEPTORS
  void InitRecept(geometry3D *GDomain,std::string nFile,int nrecep,int nt);

  // Adjoint Source

  void InitAdj(geometry3D *GDomain,std::string nFile,int nrecep,int nt);

 // FINALIZE RECEPTORS
  void EndRecept();


  // ADD VALUES RECEPTORS
    void GetRecept(int ktime);

  void printfile(Dfloat * Var,char *nfile, int ktime);

  void loadfile(Dfloat * Var,char *nfile, int ktime);

  void file(char* NameVar,int time, int io);

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

  // PRINTS INFO ABOUT RECEPTORS AND SOURCES

  void PrintInf();

  // BOUNDARIES FOR RETROPROPAGATION

  void boundX(Dfloat *var,int side,int inout,int time,int indx);

  void boundY(Dfloat *var,int side,int inout,int time,int indx);

  void boundZ(Dfloat *var,int side,int inout,int time,int indx);

  void SaveBoundaries_V(int time);

  void SaveBoundaries_S(int time);

  void LoadBoundaries_V(int time);

  void LoadBoundaries_S(int time);


  void WriteBoundaries(int time,int nt);
  void ReadBoundaries(int time,int nt);

  


  

  





  





};

#endif

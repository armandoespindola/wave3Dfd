#ifndef __PML__
#define __PML__

#include "definitions.hpp"


class DPML {

protected:
	Dfloat K_x;
	Dfloat *d_x,*d_x_half,*alpha_x;
	Dfloat *alpha_x_half,*a_x;
	Dfloat *a_x_half;
	Dfloat *b_x,*b_x_half;
	Dfloat d0_x;
	int NXT,n_pml,n_x;
	Dfloat NPOWER,ALPHA_MAX_PML,cp,Rcoef; 
	Dfloat thickness_PML_x;
  	Dfloat X_OrgL,X_OrgR;
  	Dfloat DELTAX,val,val_NORM;
public:
	// NX NODES	WITHOUT PML
	// PML_NX	NODES WITH PML
	// T_PML 	THICKNESS PML
	// Delta 	DELTA
	// DPML 	NODES PML
	// dt 	DELTA TIME
	// f0 	SOURCE FREQUENCY
	// bool	MIN_PML
	// bool MAX_PML
	DPML(int NX, int PML_NX, Dfloat T_PML, Dfloat Delta, int IPML,Dfloat dt, \
		Dfloat f0, bool, bool);
	~DPML();

	Dfloat B(int i);

	Dfloat B_HALF(int i);

	Dfloat A(int i);

	Dfloat A_HALF(int i);










};






#endif
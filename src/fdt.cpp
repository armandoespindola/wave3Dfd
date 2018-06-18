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

#include "fdt.hpp"

DFT::DFT(SDM *in_sdm,std::string nFile,int infreq){

  nf = infreq;
  sdm = in_sdm;

  freq = new Dfloat[nf];



  // READING FREQUENCIES FILE
  std::ifstream R;
  std::string line;

  R.open(nFile.c_str());


  std::getline(R,line);
  
  for (int i=0; i<nf; ++i){

    std::getline(R,line);

    freq[i] = std::stof(line);

  }

  R.close();

  Fvx = new Dfloat*[nf];
  Fvy = new Dfloat*[nf];
  Fvz = new Dfloat*[nf];
  iFvx = new Dfloat*[nf];
  iFvy = new Dfloat*[nf];
  iFvz = new Dfloat*[nf];

  
  for (int i=0;i<nf;++i){


    Fvx[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    Fvy[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    Fvz[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFvx[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFvy[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFvz[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    
  }

  InitVar(ZERO);
  

}


DFT::~DFT(){

  for (int i=0;i<nf;++i){

    delete [] Fvx[i];
    delete [] Fvy[i];
    delete [] Fvz[i];
    delete [] iFvx[i];
    delete [] iFvy[i];
    delete [] iFvz[i];
  }

  delete [] Fvx;
  delete [] Fvy;
  delete [] Fvz;
  delete [] iFvx;
  delete [] iFvy;
  delete [] iFvz;

}


void DFT::InitVar(Dfloat f){

  for (int l =0;l<nf;++l){
    for (int i=0; i<sdm->SDMGeom->HALO_Node(); ++i){

      Fvx[l][i] = f;
      Fvy[l][i] = f;
      Fvz[l][i] = f;

      iFvx[l][i] = f;
      iFvy[l][i] = f;
      iFvz[l][i] = f;

    }
  }

}

void DFT::FD(Dfloat dt,int ktime){

  Dfloat a1,a2;

  for (int l =0;l<nf;++l){

    a1 = cos( 2 * pi * freq[l] * dt * ktime) * dt;
    a2 = sin( 2 * pi * freq[l] * dt * ktime) * dt;



#pragma omp parallel for num_threads(sdm->N_omp) \
  firstprivate(l,a1,a2)
    for (int k=0;k<sdm->SNodeZ();k++){
      for (int j=0;j<sdm->SNodeY();j++){
	for (int i=0;i<sdm->SNodeX();i++){

	  Fvx[l][sdm->IJK(i,j,k)] += sdm->vx[sdm->IJK(i,j,k)] * a1;
	  iFvx[l][sdm->IJK(i,j,k)] += sdm->vx[sdm->IJK(i,j,k)] * a2;

	  Fvy[l][sdm->IJK(i,j,k)] += sdm->vy[sdm->IJK(i,j,k)] * a1;
	  iFvy[l][sdm->IJK(i,j,k)] += sdm->vy[sdm->IJK(i,j,k)] * a2;

	  Fvz[l][sdm->IJK(i,j,k)] += sdm->vz[sdm->IJK(i,j,k)] * a1;
	  iFvz[l][sdm->IJK(i,j,k)] += sdm->vz[sdm->IJK(i,j,k)] * a2;

	  
	  
	}
      }
    }

    
  }
}

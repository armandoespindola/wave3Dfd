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

  Fux = new Dfloat*[nf];
  Fuy = new Dfloat*[nf];
  Fuz = new Dfloat*[nf];
  iFux = new Dfloat*[nf];
  iFuy = new Dfloat*[nf];
  iFuz = new Dfloat*[nf];

  
  for (int i=0;i<nf;++i){
    Fux[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    Fuy[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    Fuz[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFux[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFuy[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
    iFuz[i] = new Dfloat[sdm->SDMGeom->HALO_Node()];
  }

  InitVar(ZERO);
  

}


DFT::~DFT(){

  for (int i=0;i<nf;++i){
    delete [] Fux[i];
    delete [] Fuy[i];
    delete [] Fuz[i];
    delete [] iFux[i];
    delete [] iFuy[i];
    delete [] iFuz[i]; 
  }


  delete [] Fux;
  delete [] Fuy;
  delete [] Fuz;
  delete [] iFux;
  delete [] iFuy;
  delete [] iFuz;

}


void DFT::InitVar(Dfloat f){

  for (int l =0;l<nf;++l){
    for (int i=0; i<sdm->SDMGeom->HALO_Node(); ++i){

      Fux[l][i] = f;
      Fuy[l][i] = f;
      Fuz[l][i] = f;

      iFux[l][i] = f;
      iFuy[l][i] = f;
      iFuz[l][i] = f;

    }
  }

}

void DFT::FD(Dfloat dt,int ktime){

  double a1,a2;

  for (int l =0;l<nf;++l){

    a1 = cos( 2 * pi * freq[l] * dt * ktime) * dt;
    a2 = sin( 2 * pi * freq[l] * dt * ktime) * dt;



#pragma omp parallel for num_threads(sdm->N_omp) \
  firstprivate(l,a1,a2)
    for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){

	  Fux[l][sdm->IJK(i,j,k)] += sdm->ux[sdm->IJK(i,j,k)] * a1;
	  iFux[l][sdm->IJK(i,j,k)] += sdm->ux[sdm->IJK(i,j,k)] * a2;

	  Fuy[l][sdm->IJK(i,j,k)] += sdm->uy[sdm->IJK(i,j,k)] * a1;
	  iFuy[l][sdm->IJK(i,j,k)] += sdm->uy[sdm->IJK(i,j,k)] * a2;

	  Fuz[l][sdm->IJK(i,j,k)] += sdm->uz[sdm->IJK(i,j,k)] * a1;
	  iFuz[l][sdm->IJK(i,j,k)] += sdm->uz[sdm->IJK(i,j,k)] * a2;

	  
	  
	}
      }
    }

    
  }
}

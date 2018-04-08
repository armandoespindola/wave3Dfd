/*
                  ###  Wave3Dfd (Time-Kernels) 1.0 ####

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

#include "definitions.hpp"
#include "geometry3D.hpp"
#include "model.hpp"
#include "show.hpp"
#include "sdm.hpp"
#include "mpi_trans.hpp"
#include "parameters.hpp"
#include "kernels.hpp"
#include "fdt.hpp"


int main (int argc, char* argv[]) {


  // READ PARAMETERS FILE 

  PAR par(argc,argv,"-nFile");


  // ####### PARAMETERS #######
  // ##########################


  // PARAMETERS SUBDOMAINS
 
  const int Sx = std::stoi(par.ParamReturn("-sx"));
  const int Sy = std::stoi(par.ParamReturn("-sy"));
  const int Sz = std::stoi(par.ParamReturn("-sz"));


  // NUMBER OMP THREADS

  const int N_omp = std::stoi(par.ParamReturn("-nomp"));


  // INITIAL LIMIT DOMAIN

  const VecF Ilim = {std::stof(par.ParamReturn("-xi")), \
                     std::stof(par.ParamReturn("-yi")), \
                     std::stof(par.ParamReturn("-zi"))};          


  // END LIMIT DOMAIN

  const VecF Flim = {std::stof(par.ParamReturn("-xf")), \
                     std::stof(par.ParamReturn("-yf")), \
                     std::stof(par.ParamReturn("-zf"))};          

 // NUMBER OF ELEMENTS

  const VecI Nelem = {std::stoi(par.ParamReturn("-nx"))-1, \
                          std::stoi(par.ParamReturn("-ny"))-1, \
                          std::stoi(par.ParamReturn("-nz"))-1};


  // TIME STEP

  const Dfloat dt = std::stof(par.ParamReturn("-dt"));


  // NUMBER OF SOURCES

  const int nsource = std::stoi(par.ParamReturn("-ns"));

  // SOURCES FILE

  const std::string sourceFile = par.ParamReturn("-sf");

  // NUMBER OF RECEPTORS

  const int nrecep = std::stoi(par.ParamReturn("-nr"));

  // RECEPTORS FILE

  const std::string recepFile = par.ParamReturn("-rf");

  // SOURCE FREQUENCY

  const Dfloat f0 = std::stof(par.ParamReturn("-f0"));

  // SOURCE TIME FUNCTION

  const int s_type = std::stoi(par.ParamReturn("-stype"));


  // SIMULATION TIME

  const Dfloat t = std::stof(par.ParamReturn("-t"));

  // SNAPSOTH

  const int t_snap = std::stoi(par.ParamReturn("-t_snap"));
  const int snap = std::stoi(par.ParamReturn("-snap"));

  const int t_snapR = std::stoi(par.ParamReturn("-t_snapR"));
  const int snapR = std::stoi(par.ParamReturn("-snapR"));

  // ADJOINT PROPAGATION
  const int ADJ_P = std::stoi(par.ParamReturn("-adj"));

  
  

  int  nt = round(t / dt);                                     // TIME MARCHING STEP
  VecI SubN = {Sx,Sy,Sz};                                     // NUMBER OF SUBDOMAINS
  double time1,time2,time;


  SDM *sdm,*RTP,*ADJ;                          // Pointer SubDomains
  Dfloat *SubMu,*SubRho,*SubLamb;               // Model Subdomains
  int N_mpi = SubN.x * SubN.y * SubN.z;         // Number MPI processors
  geometry3D *Gdomain;                          // Domain
  Show show;                                    // Print Class
  MPI_DATA *HALO,*HALO_ADJ;                     // Subdoamin HALO Transfer Class
  MODEL *model;                                 // Model
  DFT *uw_sdm, *uw_adj;
 


// General Domain
  Gdomain = new geometry3D(Ilim,Flim,Nelem); 

// Local nodes domain
  VecI GNod = {Gdomain->L_NodeX(),Gdomain->L_NodeY(),Gdomain->L_NodeZ()}; 

// Dimension Subdomains Nodes
  VecI SubNodes = {Gdomain->HALO_NodeX() / SubN.x, \
		   Gdomain->HALO_NodeY() / SubN.y, \
		   Gdomain->HALO_NodeZ() / SubN.z};


// Creation of SubDomain Index
VecI subi[N_mpi];

for(int k=0; k<SubN.z; ++k){
  for(int j=0; j<SubN.y; ++j){
    for(int i=0; i<SubN.x; ++i){
      int indx = i + j * SubN.x + k * SubN.y * SubN.x; 
      subi[indx] = {i,j,k};
    }
  }
}


// LIMITS OF DOMAIN WITH PML

VecF GI = {Gdomain->CoorX(0),Gdomain->CoorY(0),Gdomain->CoorZ(0)};
VecF GF = {Gdomain->CoorX(Gdomain->HALO_NodeX()-1),\
    Gdomain->CoorY(Gdomain->HALO_NodeY()-1),\
    Gdomain->CoorZ(Gdomain->HALO_NodeZ()-1)};

    //printf("%f\t%f\t%f\n",GF.x,GF.y,GF.z);



 // LIMITS OF SUBDOAMINS
 
 VecF SGI[N_mpi],SGF[N_mpi];

 for (int nsub=0; nsub<N_mpi; ++nsub){
  SGI[nsub] = {Gdomain->CoorX(subi[nsub].x * SubNodes.x), \
    Gdomain->CoorY(subi[nsub].y * SubNodes.y), \
    Gdomain->CoorZ(subi[nsub].z * SubNodes.z)};

  SGF[nsub] = {Gdomain->Dx() * (SubNodes.x - 1) + SGI[nsub].x, \
    Gdomain->Dy() * (SubNodes.y - 1) + SGI[nsub].y, \
    Gdomain->Dz() * (SubNodes.z - 1) + SGI[nsub].z};  
    //printf("%f\t%f\t%f\t",SGI[nsub].x,SGI[nsub].y,SGI[nsub].z);
    //printf("%f\t%f\t%f\n",SGF[nsub].x,SGF[nsub].y,SGF[nsub].z);
 }     

MPI_Status status;
int rank,total_proc;

MPI::Init(argc, argv);

total_proc = MPI::COMM_WORLD.Get_size();

rank = MPI::COMM_WORLD.Get_rank();


if (total_proc != N_mpi) {
    printf("Number of MPI_PROC: %d is not equal to Number Subdomains: %d\n\n",total_proc,N_mpi);
    return 0;
}


if (rank == 0) {

  std::cout<<"Parameters Subdomain and OpenMP threads"<<std::endl;
  std::cout<<"Sx : "<<Sx<<" "<<"Sy : "<<Sy<<" "<<"Sz : "<<Sx<<std::endl;
  std::cout<<"N_omp : "<<N_omp<<std::endl;
  std::cout<<"DX "<<Gdomain->Dx()<<std::endl;
  std::cout<<"DY "<<Gdomain->Dy()<<std::endl;
  std::cout<<"DZ "<<Gdomain->Dz()<<std::endl;
  
  std::cout<<"HALO_NodX "<<Gdomain->HALO_NodeX()<<std::endl;
  std::cout<<"HALO_NodY "<<Gdomain->HALO_NodeY()<<std::endl;
  std::cout<<"HALO_NodZ "<<Gdomain->HALO_NodeZ()<<std::endl;

  std::cout<<"Nodes_Subdomain X "<<SubNodes.x<<std::endl;
  std::cout<<"Nodes_Subdomain Y "<<SubNodes.y<<std::endl;
  std::cout<<"Nodes_Subdomain Z "<<SubNodes.z<<std::endl;
  std::cout<<"Time steps: "<<nt<<std::endl;
  std::cout<<"Source File: "<<sourceFile<<std::endl;
 

  

// Read Model Domain
  std::string FileVP = par.ParamReturn("-VP_F"); // "../src/example/VP.bin";
  std::string FileVS = par.ParamReturn("-VS_F"); // "../src/example/VS.bin";
  std::string FileR =  par.ParamReturn("-RHO_F"); //"../src/example/RHO.bin";

   std::cout<<"File VP: "<<FileVP<<std::endl;
   std::cout<<"File VS: "<<FileVS<<std::endl;
   std::cout<<"File RHO: "<<FileR<<std::endl;

  model = new MODEL(FileVP.c_str(),FileVS.c_str(),FileR.c_str(), \
        GNod,SubNodes);

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
  

// OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);
  

 // SubDomain Model Parameters 
  SubMu = new Dfloat[sdm->SNodeT()];
  SubRho = new Dfloat[sdm->SNodeT()];
  SubLamb = new Dfloat[sdm->SNodeT()];


  for (int i=total_proc-1; i>=0; --i){

  model->SubModel(subi[i],SubRho,SubMu,SubLamb);

  if ( i > 0) {

  MPI::COMM_WORLD.Send(SubRho,sdm->SNodeT(),MY_MPI_Dfloat,i,0);
  MPI::COMM_WORLD.Send(SubMu,sdm->SNodeT(),MY_MPI_Dfloat,i,0);
  MPI::COMM_WORLD.Send(SubLamb,sdm->SNodeT(),MY_MPI_Dfloat,i,0);

  }

  }

  //model->SubModel(subi[0],SubRho,SubMu,SubLamb);

 } 

if (rank > 0) {

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
  

  // OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);

 // SubDomain Model Parameters 
  SubMu = new Dfloat[sdm->SNodeT()];
  SubRho = new Dfloat[sdm->SNodeT()];
  SubLamb = new Dfloat[sdm->SNodeT()];

  MPI::COMM_WORLD.Recv(SubRho,sdm->SNodeT(),MY_MPI_Dfloat,0,0);
  MPI::COMM_WORLD.Recv(SubMu,sdm->SNodeT(),MY_MPI_Dfloat,0,0);
  MPI::COMM_WORLD.Recv(SubLamb,sdm->SNodeT(),MY_MPI_Dfloat,0,0);

 }

MPI_Barrier(MPI_COMM_WORLD); 


 // AUXUILARY VARIABLES
  
   char name[100];
   int NXT = Gdomain->HALO_NodeX();
   int NYT = Gdomain->HALO_NodeY();
   int NZT = Gdomain->HALO_NodeZ();




  //#########################
  // SUBDOMAIM PROGRAM 
  //#########################


// MODEL INITIALIZATION
 
  sdm->ModelRead(SubRho,"RHO");
  sdm->ModelRead(SubMu,"MU");
  sdm->ModelRead(SubLamb,"LAMB");

  // SOURCE INITIALIZATION
  sdm->InitSource(Gdomain,sourceFile,nsource);
  sdm->InitRecept(Gdomain,recepFile,nrecep,nt);
  
  int a=sdm->CFL();
  sdm->InitVar(ZERO);
  HALO = new MPI_DATA(sdm);

  time = 0.0;

  int kk = t_snap;


  // FREQUENCY DOMAIN WAVEFIELD
  //  uw_sdm = new DFT(sdm,freqFile,nfreq);

  // #############################
  // #############################
  // LOOP TIME FORWARD PROPAGATION
  // #############################
  // #############################
  
  for (int k = 0; k<nt; ++k){


    time1 = MPI::Wtime();

     if (rank == 0){
      printf("Time step : %d of %d rank %d tproc %d\n",k,nt,rank,total_proc);
    }
     // Source #########
     sdm->AddSource(k,s_type);
     sdm->GetRecept();

     if (ADJ_P)
     sdm->SaveBoundaries_S(k);

    MPI_Barrier(MPI_COMM_WORLD);
    
    HALO->TRANSFER("SXX");
    HALO->TRANSFER("SYY");
    HALO->TRANSFER("SZZ");
    HALO->TRANSFER("SXY");
    HALO->TRANSFER("SXZ");
    HALO->TRANSFER("SYZ");

    MPI_Barrier(MPI_COMM_WORLD); 

    if (ADJ_P)
    sdm->SaveBoundaries_V(k);
	
    sdm->FDVX();
    sdm->FDVY();
    sdm->FDVZ();

     HALO->TRANSFER("VX");
     HALO->TRANSFER("VY");
     HALO->TRANSFER("VZ");

    MPI_Barrier(MPI_COMM_WORLD);

    sdm->FDSII();
    sdm->FDSXY();
    sdm->FDSXZ();
    sdm->FDSYZ();

    
    time2 = MPI::Wtime();

    time += (time2-time1)/ (Dfloat) nt;




    // SAVE THE LAST BOUNDARY FOR RETROPROPAGATION
    if (ADJ_P) {
      if (k == nt-1){

	sdm->file("VZ",k,0);
	sdm->file("VY",k,0);
	sdm->file("VX",k,0);
	sdm->file("UZ",k,0);
	sdm->file("UY",k,0);
	sdm->file("UX",k,0);
	sdm->file("SXX",k,0);
	sdm->file("SYY",k,0);
	sdm->file("SZZ",k,0);
	sdm->file("SXY",k,0);
	sdm->file("SXZ",k,0);
	sdm->file("SYZ",k,0);
      }
    }
    

    
    if (snap){ 

      if (kk == k){

	  sprintf(name,"DATA/VX-%d.bin",k);
	  HALO->MergePrint(sdm->vx,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/VY-%d.bin",k);
	  HALO->MergePrint(sdm->vy,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/VZ-%d.bin",k);
	  HALO->MergePrint(sdm->vz,NXT,NYT,NZT,subi,rank,name);
      
	  kk += t_snap;
      }
    }
    
    
  } // TIME FORWARD PROPAGATION


  if (rank == 0){
    printf("END\t TIME(STEP): %f\t HOLE TIME: %f\n",time,time*nt);
  }

  sdm->EndRecept();
  delete sdm;  


  if (ADJ_P) {

  delete  HALO;
  // ###########################
  // ###########################
  // RETROPROPAGATION AND ADJOINT
  // ###########################
  // ###########################


  if (rank == 0) {

    // RETRO-PROPAGATION
    RTP = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,1);
    RTP->set_omp(N_omp);

    // ADJOINT-PROPAGATION

    ADJ = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
    ADJ->set_omp(N_omp);    



    for (int i=total_proc-1; i>=0; --i){

      model->SubModel(subi[i],SubRho,SubMu,SubLamb);

      if ( i > 0) {

	MPI::COMM_WORLD.Send(SubRho,RTP->SNodeT(),MY_MPI_Dfloat,i,0);
	MPI::COMM_WORLD.Send(SubMu,RTP->SNodeT(),MY_MPI_Dfloat,i,0);
	MPI::COMM_WORLD.Send(SubLamb,RTP->SNodeT(),MY_MPI_Dfloat,i,0);
      }
    }
    
    
  } else if (rank > 0) {


    // RETRO-PROPAGATION
    RTP = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,1);
    RTP->set_omp(N_omp);

    
    // ADJOINT-PROPAGATION    
    ADJ = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
    ADJ->set_omp(N_omp);   

    MPI::COMM_WORLD.Recv(SubRho,RTP->SNodeT(),MY_MPI_Dfloat,0,0);
    MPI::COMM_WORLD.Recv(SubMu,RTP->SNodeT(),MY_MPI_Dfloat,0,0);
    MPI::COMM_WORLD.Recv(SubLamb,RTP->SNodeT(),MY_MPI_Dfloat,0,0);
    
  }


  MPI_Barrier(MPI_COMM_WORLD); 

  

  //#########################
  // SUBDOMAIM PROGRAM 
  //#########################

  // KERNELS

  KERNEL *KS;

  KS = new KERNEL(RTP,ADJ);


// MODEL INITIALIZATION

  // RETRO-PROPAGATION
  RTP->ModelRead(SubRho,"RHO");
  RTP->ModelRead(SubMu,"MU");
  RTP->ModelRead(SubLamb,"LAMB");

  // ADJOINT-PROPAGATION
  ADJ->ModelRead(SubRho,"RHO");
  ADJ->ModelRead(SubMu,"MU");
  ADJ->ModelRead(SubLamb,"LAMB");

  // SOURCE INITIALIZATION

  // RETRO-PROPAGATION
  RTP->InitSource(Gdomain,sourceFile,nsource);
  RTP->InitVar(ZERO);

  
  // ADJOINT-PROPAGATION
  ADJ->InitAdj(Gdomain,recepFile,nrecep,nt);
  ADJ->InitVar(ZERO);

  
  time = 0.0;
  
  int Rkk = nt-t_snapR;


  // LOAD LAST BOUNDARY
  
  RTP->file("VZ",nt-1,1);
  RTP->file("VY",nt-1,1);
  RTP->file("VX",nt-1,1);
  RTP->file("UZ",nt-1,1);
  RTP->file("UY",nt-1,1);
  RTP->file("UX",nt-1,1);
  RTP->file("SXX",nt-1,1);
  RTP->file("SYY",nt-1,1);
  RTP->file("SZZ",nt-1,1);
  RTP->file("SXY",nt-1,1);
  RTP->file("SXZ",nt-1,1);
  RTP->file("SYZ",nt-1,1);
  char cleanfile[100];
  sprintf(cleanfile,"rm *%d-*%d.bin",rank,nt-1);
  system (cleanfile); 

  //RETRO-PROPAGATION TRANSFER DATA MPI
  HALO = new MPI_DATA(RTP);

  //ADJOINT TRANSFER DATA MPI
  HALO_ADJ = new MPI_DATA(ADJ);

  // #############################
  // #############################
  // LOOP TIME RETRO AND ADJOINT PROPAGATION
  // #############################
  // #############################
  
  for (int k = nt-1; k>0; --k){


    time1 = MPI::Wtime();
    if (rank == 0){
      printf("Time step (ADJ) : %d of %d rank %d tproc %d\n",k,nt,rank,total_proc);
    }


    //####################
    //###################
    // PUT ADJOINT HERE
    // THAT ENSURES THAT WE ARE DO IN THE CORRELATION AT THE SAME TIME
     
     
    MPI_Barrier(MPI_COMM_WORLD);
    
    HALO_ADJ->TRANSFER("SXX");
    HALO_ADJ->TRANSFER("SYY");
    HALO_ADJ->TRANSFER("SZZ");
    HALO_ADJ->TRANSFER("SXY");
    HALO_ADJ->TRANSFER("SXZ");
    HALO_ADJ->TRANSFER("SYZ");

    MPI_Barrier(MPI_COMM_WORLD); 
	
    ADJ->FDVX();
    ADJ->FDVY();
    ADJ->FDVZ();

     //  ADJOINT SOURCE
     ADJ->AddSourceAdj(k);
    
     HALO_ADJ->TRANSFER("VX");
     HALO_ADJ->TRANSFER("VY");
     HALO_ADJ->TRANSFER("VZ");

    MPI_Barrier(MPI_COMM_WORLD);

    ADJ->FDSII();
    ADJ->FDSXY();
    ADJ->FDSXZ();
    ADJ->FDSYZ();

    
    if (snapR){ 

      if (Rkk == k){

	  sprintf(name,"DATA/AVX-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vx,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/AVY-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vy,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/AVZ-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vz,NXT,NYT,NZT,subi,rank,name);
      

      }
    }
    
    // ^^^^^^^^^^^^^^^^^^^
    // ADJOINT PROPAGATION 
    // ###################


    //################ KERNELS ###################

    KS->CALC();


    if (snapR){ 
      if (Rkk == k){

	char name[100];
	  
	  int NXT = Gdomain->HALO_NodeX();
	  int NYT = Gdomain->HALO_NodeY();
	  int NZT = Gdomain->HALO_NodeZ();
	
	 sprintf(name,"DATA/RVX-%d.bin",k);
	  HALO->MergePrint(RTP->vx,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/RVY-%d.bin",k);
	  HALO->MergePrint(RTP->vy,NXT,NYT,NZT,subi,rank,name);

	  sprintf(name,"DATA/RVZ-%d.bin",k);
	  HALO->MergePrint(RTP->vz,NXT,NYT,NZT,subi,rank,name);
	  
      Rkk -= t_snapR;
      }
    }
    
    //printf("%f\n",source);
    RTP->FDSII();
    RTP->FDSXY();
    RTP->FDSXZ();
    RTP->FDSYZ();
    
    RTP->LoadBoundaries_S(k);

    MPI_Barrier(MPI_COMM_WORLD);
    
    HALO->TRANSFER("SXX");
    HALO->TRANSFER("SYY");
    HALO->TRANSFER("SZZ");
    HALO->TRANSFER("SXY");
    HALO->TRANSFER("SXZ");
    HALO->TRANSFER("SYZ");

    MPI_Barrier(MPI_COMM_WORLD); 
    
    RTP->FDVX();
    RTP->FDVY();
    RTP->FDVZ();
    
    RTP->LoadBoundaries_V(k);

     HALO->TRANSFER("VX");
     HALO->TRANSFER("VY");
     HALO->TRANSFER("VZ");



     // Source #########
     RTP->AddSource(k,s_type);


    MPI_Barrier(MPI_COMM_WORLD); 

    time2 = MPI::Wtime();

    time += (time2-time1)/(Dfloat)nt;    
    

  }


  KS->KERNELS();


  int NXT = Gdomain->HALO_NodeX();
  int NYT = Gdomain->HALO_NodeY();
  int NZT = Gdomain->HALO_NodeZ();

  HALO->MergePrint(KS->KDEN,NXT,NYT,NZT,subi,rank,"DATA/KRHO.bin");
  HALO->MergePrint(KS->KVS,NXT,NYT,NZT,subi,rank,"DATA/KVS.bin");
  HALO->MergePrint(KS->KVP,NXT,NYT,NZT,subi,rank,"DATA/KVP.bin");


  if (rank == 0){
    printf("END\t TIME(STEP) ADJOINT: %f\t HOLE TIME: %f\n",time,time*nt);
  }

  

  if (rank == 0){
    delete model;
    system("rm temp/*");
  }
  
  delete RTP;
  delete ADJ;
  delete HALO_ADJ;

  }// Adjoint IF
  

  delete [] SubMu;
  delete [] SubRho;
  delete [] SubLamb;
  delete  HALO;
  //delete uw_sdm;

   
MPI::Finalize();
return 0;

}

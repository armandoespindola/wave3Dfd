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

#include "wave3Dfd.hpp"

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

  // SNAPSHOT

  const int t_snap = std::stoi(par.ParamReturn("-t_snap"));
  const int snap = std::stoi(par.ParamReturn("-snap"));

  const int t_snapR = std::stoi(par.ParamReturn("-t_snapR"));
  const int snapR = std::stoi(par.ParamReturn("-snapR"));

  // ADJOINT PROPAGATION
  const int ADJ_P = std::stoi(par.ParamReturn("-adj"));

  // FREQUENCY KERNELS
  const int nfreq = std::stoi(par.ParamReturn("-nfq"));

  const std::string freqFile = par.ParamReturn("-rfq");

   // SRCFILE
  const int SrcFile = std::stoi(par.ParamReturn("-srcFile"));
  
  

  int  nt = round(t / dt);                                     // TIME MARCHING STEP
  VecI SubN = {Sx,Sy,Sz};                                     // NUMBER OF SUBDOMAINS
  double time1,time2,time;


  SDM *sdm,*RTP,*ADJ;                           // Pointer SubDomains
  Dfloat *SubMod;                               // Model Subdomains
  int N_mpi = SubN.x * SubN.y * SubN.z;         // Number MPI processors
  geometry3D *Gdomain;                          // Domain
  Show show;                                    // Print Class
  MPI_DATA *HALO,*HALO_ADJ;                     // Subdoamin HALO Transfer Class
  MODEL *model;                                 // Model
  DFT *uw_sdm, *uw_adj;
  char name[100];
  int NXT,NYT, NZT;
  int cfl;


 


// General Domain
  Gdomain = new geometry3D(Ilim,Flim,Nelem); 

// Local nodes domain
  VecI GNod = {Gdomain->L_NodeX(),Gdomain->L_NodeY(),Gdomain->L_NodeZ()}; 

// Dimension Subdomains Nodes
  VecI SubNodes = {Gdomain->HALO_NodeX() / SubN.x, \
		   Gdomain->HALO_NodeY() / SubN.y, \
		   Gdomain->HALO_NodeZ() / SubN.z};

  MPI_Status status;
  int rank,total_proc;
  int ndims = 3;
  int dim[3],coords[3];
  int isperiodic[3];
  int reorder = 1;
  MPI_Comm COM3D;
  
  isperiodic[0] = 0;
  isperiodic[1] = 0;
  isperiodic[2] = 0;
  
  dim[0] = SubN.x;
  dim[1] = SubN.y;
  dim[2] = SubN.z;

  //MPI::Init(argc, argv);
  MPI_Init(&argc,&argv);

  //total_proc = MPI::COMM_WORLD.Get_size();
  MPI_Comm_size(MPI_COMM_WORLD,&total_proc);

  // Define 3D MPI DOMAIN
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dim
		 ,isperiodic, reorder, &COM3D);


  //rank = MPI::COMM_WORLD.Get_rank();
  MPI_Comm_rank(COM3D,&rank);

  MPI_Cart_coords(COM3D, rank, ndims,coords);

  
  //std::cout<<rank<<" "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
				

// Creation of SubDomain Index
 VecI subi[N_mpi];

 subi[rank] = {coords[0],coords[1],coords[2]};

 MPI_Allgather(&subi[rank],3,MPI_INT,&subi,3,MPI_INT,COM3D);


 // if (rank==5) {
 //   for(int k=0;k<N_mpi;k++){
 //     std::cout<<k<<" "<<subi[k].x<<"-"<<subi[k].y<<"-"<<subi[k].z<<std::endl;
 //   }
 // }
 
// // Creation of SubDomain Index
// VecI subi[N_mpi];

// for(int k=0; k<SubN.z; ++k){
//   for(int j=0; j<SubN.y; ++j){
//     for(int i=0; i<SubN.x; ++i){
//       int indx = i + j * SubN.x + k * SubN.y * SubN.x; 
//       subi[indx] = {i,j,k};
//     }
//   }
// }


// LIMITS OF DOMAIN WITH PML

VecF GI = {Gdomain->CoorX(0),Gdomain->CoorY(0),Gdomain->CoorZ(0)};
VecF GF = {Gdomain->CoorX(Gdomain->HALO_NodeX()-1), \
	   Gdomain->CoorY(Gdomain->HALO_NodeY()-1), \
	   Gdomain->CoorZ(Gdomain->HALO_NodeZ()-1)};
//printf("%f\t%f\t%f\n",GF.x,GF.y,GF.z);



 // LIMITS OF SUBDMAINS
 
 VecF SGI[N_mpi],SGF[N_mpi];


 SGI[rank] = {Gdomain->CoorX(subi[rank].x * SubNodes.x), \
	      Gdomain->CoorY(subi[rank].y * SubNodes.y), \
	      Gdomain->CoorZ(subi[rank].z * SubNodes.z)};


 SGF[rank] = {Gdomain->Dx() * (SubNodes.x - 1) + SGI[rank].x, \
	      Gdomain->Dy() * (SubNodes.y - 1) + SGI[rank].y, \
	      Gdomain->Dz() * (SubNodes.z - 1) + SGI[rank].z};  

   //printf("%f\t%f\t%f\t",SGI[nsub].x,SGI[nsub].y,SGI[nsub].z);
   //printf("%f\t%f\t%f\n",SGF[nsub].x,SGF[nsub].y,SGF[nsub].z); 



if (total_proc != N_mpi) {
    if (rank==0){
    printf("Number of MPI_PROC: %d is not equal to Number Subdomains: %d\n",total_proc,N_mpi);
     printf("EDIT  parameter.par File\n");
  }
    MPI_Finalize();
    return 0;
}

time1 = MPI_Wtime();

if (rank == 0) {

  std::cout<<"Parameters Subdomain and OpenMP threads"<<std::endl;
  std::cout<<"Sx : "<<Sx<<" "<<"Sy : "<<Sy<<" "<<"Sz : "<<Sz<<std::endl;
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

  cfl = model->CFL(dt,Gdomain->Dx(),Gdomain->Dy(),Gdomain->Dz());

 }

MPI_Bcast(&cfl, 1, MPI_INT, 0, COM3D);

// CFL Condition
 
 if (cfl != 1){
   if (rank == 0){
     if (cfl == 0)
       printf("CFL NOT SATISFIED\n");
     if (cfl == -1)
       printf("ERROR IN MODEL PARAMETERS\n");
     delete sdm;
     delete model;
     
   }
   delete Gdomain;
   MPI_Finalize();
   return 0;
 }

 


 if (rank == 0) { 

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
  

// OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);
  

 // SubDomain Model Parameters 
  SubMod = new Dfloat[sdm->SNodeT() * 3];

  for (int i=total_proc-1; i>=0; --i){

// SubRho SubMu SubLamb

  model->SubModel(subi[i],SubMod,SubMod + sdm->SNodeT() ,SubMod + 2 * sdm->SNodeT());

  if ( i > 0) {

//  MPI::COMM_WORLD.Send(SubRho,sdm->SNodeT(),MY_MPI_Dfloat,i,0);	
//  MPI::COMM_WORLD.Send(SubMu,sdm->SNodeT(),MY_MPI_Dfloat,i,0);
//  MPI::COMM_WORLD.Send(SubLamb,sdm->SNodeT(),MY_MPI_Dfloat,i,0);

  MPI_Send(SubMod,sdm->SNodeT() * 3,MY_MPI_Dfloat,i,0,COM3D);	


  }

  }

  //model->SubModel(subi[0],SubRho,SubMu,SubLamb);

 }

if (rank > 0) {

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
  

  // OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);

 // SubDomain Model Parameters 
  SubMod = new Dfloat[sdm->SNodeT() * 3];

//  MPI::COMM_WORLD.Recv(SubRho,sdm->SNodeT(),MY_MPI_Dfloat,0,0);
//  MPI::COMM_WORLD.Recv(SubMu,sdm->SNodeT(),MY_MPI_Dfloat,0,0);
//  MPI::COMM_WORLD.Recv(SubLamb,sdm->SNodeT(),MY_MPI_Dfloat,0,0);

  MPI_Recv(SubMod,sdm->SNodeT() * 3,MY_MPI_Dfloat,0,0,COM3D,&status);

 }


time2 = MPI_Wtime();

if (rank == 0) {
      printf("Time broadcasting model to MPI cores : %f\n",time2-time1);
    }

MPI_Barrier(COM3D); 


 // AUXUILARY VARIABLES
   NXT = Gdomain->HALO_NodeX();
   NYT = Gdomain->HALO_NodeY();
   NZT = Gdomain->HALO_NodeZ();




  //#########################
  // SUBDOMAIM PROGRAM 
  //#########################


// MODEL INITIALIZATION
 
  sdm->ModelRead(SubMod,"RHO");
  sdm->ModelRead(SubMod + sdm->SNodeT(),"MU");
  sdm->ModelRead(SubMod + 2 * sdm->SNodeT(),"LAMB");

  // SOURCE INITIALIZATION
  sdm->InitSource(Gdomain,sourceFile,nsource,SrcFile,nt);
  sdm->InitRecept(Gdomain,recepFile,nrecep,nt);

  // PRINT INFORMATION
  if (rank==0){
    sdm->PrintInf();
    std::cout<<"\n"<<std::endl;
    std::cout<<"###############################"<<std::endl;
    std::cout<<"###############################"<<std::endl;
    std::cout<<"      SIMULATION STARTS\n"<<std::endl;
    
  }
  sdm->InitVar(ZERO);
  HALO = new MPI_DATA(sdm);

  //std::cout<<"HOLA\n"<<std::endl;
  time = 0.0;

  int kk = t_snap;
  int tinf = t_snap;
   

  // #############################
  // #############################
  // LOOP TIME FORWARD PROPAGATION
  // #############################
  // #############################
  
  for (int k = 0; k<nt; ++k){


    time1 = MPI_Wtime();

     // Source #########
     sdm->AddSource(k,s_type);

    
    MPI_Barrier(COM3D);
    
    // TRANSFER STRESESS
    HALO->TRANSFER(2,COM3D);
  

    // SAVE BOUNDARIES TO RECONSTRUCT WAVEFIELD

    if (ADJ_P){
     sdm->SaveBoundaries_S(k);
     sdm->SaveBoundaries_V(k);
     sdm->WriteBoundaries(k+1,nt);
    }
    MPI_Barrier(COM3D); 
	
    sdm->FDVX();
    sdm->FDVY();
    sdm->FDVZ();

    // TRANSFER VELOCITIES
    HALO->TRANSFER(1,COM3D);

    sdm->GetRecept(k);

    // if (ADJ_P){
    // sdm->SaveBoundaries_V(k);
    // sdm->WriteBoundaries(k+1,nt);
    // }
    
    MPI_Barrier(COM3D);

    
    sdm->FDSII();
    sdm->FDSXY();
    sdm->FDSXZ();
    sdm->FDSYZ();

    
    time2 = MPI_Wtime();

    if ((rank == 0) && (tinf == k)){
      //printf("Forward Propagation ->  Time step : %d of %d - MPI_Time %f - Nproc %d\n",k,nt,(time2-time1)*t_snap,total_proc);
      std::cout<<"Forward Propagation ->  Time step: "<< \
	k<<" of "<<nt<<" - MPI_Time "<<(time2-time1)*t_snap<< \
	" - Nproc "<<total_proc<<std::endl;
      tinf += t_snap;
    }

    time += (time2-time1)/ (Dfloat) nt;
    
    if (snap){ 

      if (kk == k){

	  sprintf(name,"DATA/VX-%d.bin",k);
	  HALO->MergePrint(sdm->vx,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/VY-%d.bin",k);
	  HALO->MergePrint(sdm->vy,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/VZ-%d.bin",k);
	  HALO->MergePrint(sdm->vz,NXT,NYT,NZT,subi,rank,name,COM3D);
      
	  kk += t_snap;
      }
    }




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
    
    
  } // TIME FORWARD PROPAGATION


  
  if (rank == 0){
    printf("END\t TIME(STEP): %f\t HOLE TIME: %f\n",time,time*nt);
  }

  sdm->EndRecept();
  sdm->EndSource();
  delete sdm;
  delete HALO;



  if (ADJ_P) {



  // ############################
  // ############################
  // RETROPROPAGATION AND ADJOINT
  // ############################
  // ############################

    
    // ADJOINT-PROPAGATION

    ADJ = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,0);
    ADJ->set_omp(N_omp);

    // RETRO-PROPAGATION
    RTP = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN,1);
    RTP->set_omp(N_omp);

  

  //#########################
  // SUBDOMAIM PROGRAM 
  //#########################


  // KERNELS

  KERNEL *KS;

  KS = new KERNEL(RTP,ADJ);
  
// MODEL INITIALIZATION

  
  // ADJOINT-PROPAGATION
  ADJ->ModelRead(SubMod,"RHO");
  ADJ->ModelRead(SubMod + ADJ->SNodeT(),"MU");
  ADJ->ModelRead(SubMod + 2 * ADJ->SNodeT(),"LAMB");


  // RETRO-PROPAGATION
  RTP->ModelRead(SubMod,"RHO");
  RTP->ModelRead(SubMod + ADJ->SNodeT(),"MU");
  RTP->ModelRead(SubMod + 2 * ADJ->SNodeT(),"LAMB");

  // SOURCE INITIALIZATION
  
  // ADJOINT-PROPAGATION
  ADJ->InitAdj(Gdomain,recepFile,nrecep,nt);
  ADJ->InitVar(ZERO);

  // RETRO-PROPAGATION
  RTP->InitSource(Gdomain,sourceFile,nsource,SrcFile,nt);
  RTP->InitVar(ZERO);

  
  time = 0.0;
  
  int Rkk = nt-t_snapR;
  int tinf_adj = nt-t_snapR;



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
  

  // RETRO-PROPAGATION TRANSFER DATA MPI
  HALO = new MPI_DATA(RTP);

  HALO->TRANSFER(2,COM3D);
  HALO->TRANSFER(1,COM3D);


  // ADJOINT TRANSFER DATA MPI
  HALO_ADJ = new MPI_DATA(ADJ);

  // #############################
  // #############################
  // LOOP TIME RETRO AND ADJOINT PROPAGATION
  // #############################
  // #############################
  
  for (int k = nt-1; k>0; --k){


    time1 = MPI_Wtime();


    //####################
    //###################
    // PUT ADJOINT HERE
    // THAT ENSURES THAT WE ARE DOING THE CORRELATION AT THE SAME TIME
     
     
    MPI_Barrier(COM3D);
    

    // TRANSFER STRESESS
    HALO_ADJ->TRANSFER(2,COM3D);

    MPI_Barrier(COM3D); 
	
    ADJ->FDVX();
    ADJ->FDVY();
    ADJ->FDVZ();

     //  ADJOINT SOURCE
     ADJ->AddSourceAdj(k);
    
    // TRANSFER VELOCITIES
     HALO_ADJ->TRANSFER(1,COM3D);

    MPI_Barrier(COM3D);
    
    ADJ->FDSII();
    ADJ->FDSXY();
    ADJ->FDSXZ();
    ADJ->FDSYZ();

    
    if (snapR){ 

      if (Rkk == k){

	  sprintf(name,"DATA/AVX-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vx,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/AVY-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vy,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/AVZ-%d.bin",k);
	  HALO_ADJ->MergePrint(ADJ->vz,NXT,NYT,NZT,subi,rank,name,COM3D);

	  //Rkk -= t_snapR;

      }
    }




    // ##### KERNELS COMPUTATION ####
    KS->CALC();


    if (snapR){ 
      if (Rkk == k){

	char name[100];
	  
	  int NXT = Gdomain->HALO_NodeX();
	  int NYT = Gdomain->HALO_NodeY();
	  int NZT = Gdomain->HALO_NodeZ();
	
	 sprintf(name,"DATA/RVX-%d.bin",k);
	 HALO->MergePrint(RTP->vx,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/RVY-%d.bin",k);
	  HALO->MergePrint(RTP->vy,NXT,NYT,NZT,subi,rank,name,COM3D);

	  sprintf(name,"DATA/RVZ-%d.bin",k);
	  HALO->MergePrint(RTP->vz,NXT,NYT,NZT,subi,rank,name,COM3D);

	  
      Rkk -= t_snapR;
      }
    }
    
    //printf("%f\n",source);
    RTP->FDSII();
    RTP->FDSXY();
    RTP->FDSXZ();
    RTP->FDSYZ();

    RTP->ReadBoundaries(k+1,nt);
    RTP->LoadBoundaries_S(k);

    MPI_Barrier(COM3D);
    
    // TRANSFER STRESESS
    HALO->TRANSFER(2,COM3D);

    MPI_Barrier(COM3D); 
    
    RTP->FDVX();
    RTP->FDVY();
    RTP->FDVZ();
    
    RTP->LoadBoundaries_V(k);

    MPI_Barrier(COM3D);


    // Source #########
    RTP->AddSource(k,s_type);
     
    // TRANSFER VELOCITIES
    HALO->TRANSFER(1,COM3D);    

    MPI_Barrier(COM3D); 

    time2 = MPI_Wtime();

    if ((rank == 0) && (tinf_adj == k)){
      //printf("Adjoint Simulation -> Time step : %d of %d - MPI_Time %f - Nproc %d\n",k,nt,(time2-time1)*t_snap,total_proc);
      std::cout<<"Adjoint Propagation ->  Time step: "<< \
	k<<" of "<<nt<<" - MPI_Time "<<(time2-time1)*t_snap<< \
	" - Nproc "<<total_proc<<std::endl;
      tinf_adj -= t_snap;
    }

    time += (time2-time1)/(Dfloat)nt;    
    

  } // TIME ADJOINT



  if (rank == 0){
    printf("END\t TIME(STEP) ADJOINT: %f\t HOLE TIME: %f\n",time,time*nt);
  }


  // SAVING KERNELS
   KS->KERNELS();
   Dfloat *KER_L;	
   KER_L = new Dfloat[ADJ->SNodeT() * 5];
   int noffset = ADJ->SNodeT();

   KS->GET_K(KER_L, KER_L + noffset, KER_L + 2 * noffset);

   HALO_ADJ->KernPrint(KER_L,NXT,NYT,NZT,subi,rank,"DATA/",COM3D);

   MPI_Barrier(COM3D);
   
   delete [] KER_L;

   
    
  delete RTP;
  delete ADJ;
  delete HALO_ADJ;
  delete HALO;
  
  } // IF ADJOINT


  delete [] SubMod;

  // Clean Temporal Files

  //if (rank==0){
  //char cleanfile[] = "find ./temp -type f -delete";
  //system (cleanfile);
  //}

  MPI_Finalize();
  return 0;



}

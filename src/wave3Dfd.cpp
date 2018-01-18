#include "definitions.hpp"
#include "geometry3D.hpp"
#include "model.hpp"
#include "show.hpp"
#include "sdm.hpp"
#include "mpi_trans.hpp"
#include "parameters.hpp"


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


  // SOURCE FREQUENCY

  const Dfloat f0 = std::stof(par.ParamReturn("-f0"));


  // SIMULATION TIME

  const Dfloat t = std::stof(par.ParamReturn("-t")); 


  int  nt = (int) t / dt;                                     // TIME MARCHING STEP
  VecF SrcPosition = {11000/2.0,11000/2.0,11000/2.0};               // SOURCE POSITION
  VecI SubN = {Sx,Sy,Sz};                                     // NUMBER OF SUBDOMAINS
  Dfloat time1,time2,time;


  SDM *sdm;                                     // Pointer SubDomains
  Dfloat *SubMu,*SubRho,*SubLamb;               // Model Subdomains
  int N_mpi = SubN.x * SubN.y * SubN.z;         // Number MPI processors
  geometry3D *Gdomain;                          // Domain
  Show show;                                    // Print Class
  MPI_DATA *sHALO;                              // Subdoamin HALO Transfer Class
  MODEL *model;                                 // Model
 


// General Domain
  Gdomain = new geometry3D(Ilim,Flim,Nelem); 

 // Source position 
  VecI pos_src = Gdomain->FindNode(SrcPosition);  

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

  

// Read Model Domain
  char FileVP[] = "../src/example/VP.bin";
  char FileVS[] = "../src/example/VS.bin";
  char FileR[] = "../src/example/RHO.bin";

  model = new MODEL(FileVP,FileVS,FileR, \
        GNod,SubNodes);

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN); 

// OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);

 // SubDomain Model Parameters 
  SubMu = new Dfloat[sdm->SNodeT()];
  SubRho = new Dfloat[sdm->SNodeT()];
  SubLamb = new Dfloat[sdm->SNodeT()];


  for (int i=total_proc-1; i>=0; --i){

  model->SubModel(subi[i],SubRho,SubMu,SubLamb);

  if ( i > 0) {

  MPI::COMM_WORLD.Send(SubRho,sdm->SNodeT(),MPI_DOUBLE,i,0);
  MPI::COMM_WORLD.Send(SubMu,sdm->SNodeT(),MPI_DOUBLE,i,0);
  MPI::COMM_WORLD.Send(SubLamb,sdm->SNodeT(),MPI_DOUBLE,i,0);

  }

  }

  //model->SubModel(subi[0],SubRho,SubMu,SubLamb);

 } 

if (rank > 0) {

  sdm = new SDM(GI,GF,GNod,SGI[rank],SGF[rank],SubNodes,f0,dt,subi[rank],SubN); 

  // OMP NUMBER OF THREADS
  sdm->set_omp(N_omp);

 // SubDomain Model Parameters 
  SubMu = new Dfloat[sdm->SNodeT()];
  SubRho = new Dfloat[sdm->SNodeT()];
  SubLamb = new Dfloat[sdm->SNodeT()];

  MPI::COMM_WORLD.Recv(SubRho,sdm->SNodeT(),MPI_DOUBLE,0,0);
  MPI::COMM_WORLD.Recv(SubMu,sdm->SNodeT(),MPI_DOUBLE,0,0);
  MPI::COMM_WORLD.Recv(SubLamb,sdm->SNodeT(),MPI_DOUBLE,0,0);

 }

MPI_Barrier(MPI_COMM_WORLD); 




  //#########################
  // SUBDOMAIM PROGRAM 
  //#########################


  sdm->ModelRead(SubRho,"RHO");
  sdm->ModelRead(SubMu,"MU");
  sdm->ModelRead(SubLamb,"LAMB");
  int a=sdm->CFL();
  sdm->InitVar(ZERO);
  sHALO = new MPI_DATA(sdm);

  time = 0.0;

  int kk = 100;
  for (int k = 0; k<nt; ++k){


    time1 = MPI::Wtime();
    //printf("Time step : %d of %d rank %d tproc %d\n",k,nt,rank,total_proc);
  
    Dfloat source = 100000000.0 * sdm->source(1,k,0.25);
    //printf("%f\n",source);
    sdm->FDSII();
    sdm->FDSXY();
    sdm->FDSXZ();
    sdm->FDSYZ();

     // Source #########
    sdm->AddVal(pos_src,"SXX",source);
    sdm->AddVal(pos_src,"SYY",source);
    sdm->AddVal(pos_src,"SZZ",source);

    MPI_Barrier(MPI_COMM_WORLD);
    
    sHALO->TRANSFER("SXX");
    sHALO->TRANSFER("SYY");
    sHALO->TRANSFER("SZZ");
    sHALO->TRANSFER("SXY");
    sHALO->TRANSFER("SXZ");
    sHALO->TRANSFER("SYZ");

    MPI_Barrier(MPI_COMM_WORLD); 
    
    sdm->FDVX();
    sdm->FDVY();
    sdm->FDVZ();


     sHALO->TRANSFER("VX");
     sHALO->TRANSFER("VY");
     sHALO->TRANSFER("VZ");

    MPI_Barrier(MPI_COMM_WORLD); 

    time2 = MPI::Wtime();

    time += (time2-time1)/nt;

    
    if (kk == k){
      sdm->print("VZ",k);
      //sdm.print("VY",k);
      //sdm.print("SXX",k);
      kk += 100;
      }
    
    
  }


  if (rank == 0){
    delete model;
    printf("END\t TIME(STEP): %f\t HOLE TIME: %f\n",time,time*nt);
  }
  
  delete sdm;
  delete [] SubMu;
  delete [] SubRho;
  delete [] SubLamb;
  delete  sHALO;

   
MPI::Finalize();
return 0;

}

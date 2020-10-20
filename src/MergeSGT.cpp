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

#include "MergeSGT.hpp"




void MergeSGT(char *FileName,VecI SubN,VecI SubNodes,VecI *subi,int Nrank,VecI CPML, \
	      VecI dsg,int rank,int nt,MPI_Comm comm){
  
  char name[200];
  MPI_File fhw;
  int idx,idx2;
  Dfloat *GlobalArray = NULL;
  Dfloat *SortGlobalArray = NULL;
  int Nsub = Nrank;
  MPI_Status status;
  
  VecI dim_data = {int(SubNodes.x / dsg.x),\
		   int(SubNodes.y / dsg.y),\
		   int(SubNodes.z / dsg.z)};

  if (rank==0){
    std::cout<<"Nodes: "<<dim_data.x<<","<<dim_data.y<<","<<dim_data.z<<std::endl;
  }


  int size = dim_data.x * dim_data.y * dim_data.z * nt;


  Dfloat *ReadArray,*SortArray;

  ReadArray = new Dfloat[size];
  SortArray = new Dfloat[size];

  sprintf(name,"DATA/%s-%d.bin",FileName,rank);

  if (rank ==0){
  std::cout<< "Reading File: "<<FileName<<std::endl;
  }

  
  MPI_File_open(MPI_COMM_SELF,name,MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE,MPI_INFO_NULL,&fhw);
  MPI_File_read(fhw,ReadArray,size,MY_MPI_Dfloat,&status);
  MPI_File_close(&fhw);
  //MPI_File_delete(name, MPI_INFO_NULL );


  FOR_IJK(iz,0,dim_data.z,iy,0,dim_data.y,ix,0,dim_data.x)
    for (int ik =0;ik<nt;ik++){
      idx = ik * dim_data.x * dim_data.y * dim_data.z \
  	+ ix + iy * dim_data.x + iz * dim_data.x * dim_data.y;
      idx2 = ik + (ix + iy * dim_data.x + iz * dim_data.x * dim_data.y) * nt;
      SortArray[idx2] = ReadArray[idx];
    }
  END_FOR_IJK



  if (rank==0){
    GlobalArray = new Dfloat[size * Nrank];
    SortGlobalArray = new Dfloat[size * Nrank];
  }


  // ReadArray -> SortArray
  MPI_Gather(SortArray,size,MY_MPI_Dfloat,\
	     GlobalArray,size,MY_MPI_Dfloat,\
	     0,comm);


  MPI_Barrier(comm);
  
   if (rank==0){

    int ixg,iyg,izg;
    
    for (int isub = 0; isub<Nsub;++isub) {

      // std::cout<<subi[isub].x<<","<<subi[isub].y<<","<<subi[isub].z<<std::endl;
      FOR_IJK(iz,0,dim_data.z,iy,0,dim_data.y,ix,0,dim_data.x)
    	for (int ik =0;ik<nt;ik++){
	  

    	  ixg = ix + subi[isub].x * dim_data.x;
    	  iyg = iy + subi[isub].y * dim_data.y;
    	  izg = iz + subi[isub].z * dim_data.z;
	  
    	  idx = ik + ixg * nt + iyg * nt * dim_data.x * SubN.x + \
	    izg * nt * dim_data.x * dim_data.y * SubN.x * SubN.y;
	  
    	  idx2 = ik + (ix + iy * dim_data.x + iz * dim_data.x * dim_data.y) * nt;

    	  SortGlobalArray[idx] = GlobalArray[idx2 + size * isub];
	   

    	}

      END_FOR_IJK

    	}


    FILE *R;

    sprintf(name,"DATA/%s.bin",FileName);
    R=fopen(name,"wb");

    VecI s_cpml;

    s_cpml = {int(CPML.x / dsg.x),\
	      int(CPML.y / dsg.y),\
	      int(CPML.z / dsg.z)};

  
     for (int k=s_cpml.z;k<dim_data.z * SubN.z;k++){
       for (int j=s_cpml.y;j<dim_data.y * SubN.y - s_cpml.y;j++){
     	for (int i=s_cpml.x;i<dim_data.x * SubN.x - s_cpml.x;i++){

     	  idx = i * nt + j * nt * dim_data.x * SubN.x +\
     	    k * nt * dim_data.x * dim_data.y * SubN.x * SubN.y;
    	  fwrite(SortGlobalArray + idx,sizeof(Dfloat),nt,R);
         
     	}
       }
     }
  

    fclose(R);
    
    delete [] GlobalArray;
	       
  }
  
 delete [] ReadArray;
 delete [] SortArray;

};

int main (int argc, char* argv[]) {


  
  // READ PARAMETERS FILE 

  PAR par(argc,argv,"-nFile");


  // ####### PARAMETERS #######
  // ##########################


  // PARAMETERS SUBDOMAINS
 
  const int Sx = std::stoi(par.ParamReturn("-sx"));
  const int Sy = std::stoi(par.ParamReturn("-sy"));
  const int Sz = std::stoi(par.ParamReturn("-sz"));

  // Damping Profiles CPML 
  const VecI CPML = {std::stoi(par.ParamReturn("-pmlx")), \
                    std::stoi(par.ParamReturn("-pmly")), \
                    std::stoi(par.ParamReturn("-pmlz"))};
  

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
   // DOWNSAMPLING SPACE AXES PARAMETERS
  const VecI dsf = {std::stoi(par.ParamReturn("-dsf_x")), \
                          std::stoi(par.ParamReturn("-dsf_y")), \
                          std::stoi(par.ParamReturn("-dsf_z"))};
  
  const int t_snapR = std::stoi(par.ParamReturn("-t_snapR"));
  const int snapR = std::stoi(par.ParamReturn("-snapR"));

  // SIMULATION TYPE
  const int simulation_type = std::stoi(par.ParamReturn("-simulation_type"));


   // SRCFILE
  const int SrcFile = std::stoi(par.ParamReturn("-srcFile"));

  // ### STRAIN GREEN'S FUNCTION PARAMETERS ###
  // DOWNSAMPLING TIME AXIS PARAMETER
  const int dsk = std::stoi(par.ParamReturn("-dsk"));
  // DOWNSAMPLING SPACE AXES PARAMETERS
  const VecI dsg = {std::stoi(par.ParamReturn("-dsg_x")), \
                          std::stoi(par.ParamReturn("-dsg_y")), \
                          std::stoi(par.ParamReturn("-dsg_z"))};
  
  

  int  nt = round(t / dt);                      // TIME MARCHING STEP
  VecI SubN = {Sx,Sy,Sz};                       // NUMBER OF SUBDOMAINS
  double time1,time2,time;
  int N_mpi = SubN.x * SubN.y * SubN.z;         // Number MPI processors
  geometry3D *Gdomain;                          // Domain
  Show show;                                    // Print Class
  MPI_DATA *HALO;                               // Subdoamin HALO Transfer Class
  MODEL *model;                                 // Model
  char name[100];
  int NXT,NYT, NZT;



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
  
  
  if (total_proc != N_mpi) {
    if (rank==0){
    printf("Number of MPI_PROC: %d is not equal to Number Subdomains: %d\n",total_proc,N_mpi);
     printf("EDIT  parameter.par File\n");
  }
    MPI_Finalize();
    return 0;
}




  // General Domain
  Gdomain = new geometry3D(CPML,Ilim,Flim,Nelem); 

// Local nodes domain
  VecI GNod = {Gdomain->L_NodeX(),Gdomain->L_NodeY(),Gdomain->L_NodeZ()}; 

// Dimension Subdomains Nodes
  VecI SubNodes = {Gdomain->HALO_NodeX() / SubN.x, \
		   Gdomain->HALO_NodeY() / SubN.y, \
		   Gdomain->HALO_NodeZ() / SubN.z};
  
 				

// Creation of SubDomain Index
  VecI subi[N_mpi],subj[N_mpi];

  int xc,yc,zc;

  zc = int(rank / (SubN.x * SubN.y));
  yc = int((rank - zc * (SubN.x * SubN.y)) / SubN.x);
  xc = int((rank - zc * (SubN.x * SubN.y) - yc * SubN.x));
  
  subj[rank] = {xc,yc,zc};

 subi[rank] = {coords[0],coords[1],coords[2]};

 // MPI_Allgather(&subi[rank],3,MPI_INT,&subi,3,MPI_INT,COM3D);
 MPI_Allgather(&subj[rank],3,MPI_INT,&subj,3,MPI_INT,COM3D);

 //cd int isub = rank;
 //std::cout<<"subi "<<subi[isub].x<<","<<subi[isub].y<<","<<subi[isub].z<<","<<rank<<std::endl;
 //std::cout<<"subj "<<subj[isub].x<<","<<subj[isub].y<<","<<subj[isub].z<<","<<rank<<std::endl;


// LIMITS OF DOMAIN WITH PML

VecF GI = {Gdomain->CoorX(0),Gdomain->CoorY(0),Gdomain->CoorZ(0)};
VecF GF = {Gdomain->CoorX(Gdomain->HALO_NodeX()-1), \
	   Gdomain->CoorY(Gdomain->HALO_NodeY()-1), \
	   Gdomain->CoorZ(Gdomain->HALO_NodeZ()-1)};



 // LIMITS OF SUBDMAINS
 
 VecF SGI[N_mpi],SGF[N_mpi];


 SGI[rank] = {Gdomain->CoorX(subj[rank].x * SubNodes.x), \
	      Gdomain->CoorY(subj[rank].y * SubNodes.y), \
	      Gdomain->CoorZ(subj[rank].z * SubNodes.z)};


 SGF[rank] = {Gdomain->Dx() * (SubNodes.x - 1) + SGI[rank].x, \
	      Gdomain->Dy() * (SubNodes.y - 1) + SGI[rank].y, \
	      Gdomain->Dz() * (SubNodes.z - 1) + SGI[rank].z};



 if (rank == 0) {

  std::cout<<"### Parameters ###"<<std::endl;
  std::cout<<"Sx : "<<Sx<<" "<<"Sy : "<<Sy<<" "<<"Sz : "<<Sz<<std::endl;
  std::cout<<"N_omp : "<<N_omp<<std::endl;
  std::cout<<"DX "<<Gdomain->Dx()<<std::endl;
  std::cout<<"DY "<<Gdomain->Dy()<<std::endl;
  std::cout<<"DZ "<<Gdomain->Dz()<<std::endl;
  
  std::cout<<"NodX "<<Gdomain->L_NodeX()<<std::endl;
  std::cout<<"NodY "<<Gdomain->L_NodeY()<<std::endl;
  std::cout<<"NodZ "<<Gdomain->L_NodeZ()<<std::endl;
    
  std::cout<<"NodX (with PML) "<<Gdomain->HALO_NodeX()<<std::endl;
  std::cout<<"NodY (with PML) "<<Gdomain->HALO_NodeY()<<std::endl;
  std::cout<<"NodZ (with PML) "<<Gdomain->HALO_NodeZ()<<std::endl;

  std::cout<<"Nodes_Subdomain X "<<SubNodes.x<<std::endl;
  std::cout<<"Nodes_Subdomain Y "<<SubNodes.y<<std::endl;
  std::cout<<"Nodes_Subdomain Z "<<SubNodes.z<<std::endl;
  std::cout<<"Time steps: "<<nt<<std::endl;


 }


 // Read source name
 
 source *sourceM = new source(Gdomain,sourceFile,nsource);
 sprintf(name,"Hxx-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);


 MPI_Barrier(COM3D);


 sprintf(name,"Hxy-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);

  MPI_Barrier(COM3D);

 sprintf(name,"Hxz-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);

  MPI_Barrier(COM3D);
  
 sprintf(name,"Hyy-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);

  MPI_Barrier(COM3D);
  
 sprintf(name,"Hyz-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);

  MPI_Barrier(COM3D);
  
 sprintf(name,"Hzz-%s",sourceM->nameSource[0].c_str());
 MergeSGT(name,SubN,SubNodes,subj,N_mpi,CPML,	\
	  dsg,rank,nt,COM3D);


 MPI_Finalize();
 return 0;
  
}




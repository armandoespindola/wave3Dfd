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

#include "mpi_trans.hpp"


MPI_DATA::MPI_DATA(SDM *Isdm){
  sdm = Isdm;

  N_SN = KHALO * sdm->L_SNodeX() * sdm->L_SNodeZ();

  N_WE = KHALO * sdm->L_SNodeY() * sdm->L_SNodeZ();

  N_UpDown = KHALO * sdm->L_SNodeX() * sdm->L_SNodeY();  

    // HALO COMMUNICATION O OUT 1 IN 
  BS0 = new Dfloat[N_SN]; 
  BN0 = new Dfloat[N_SN];
  BS1 = new Dfloat[N_SN]; 
  BN1 = new Dfloat[N_SN];
  
  BW0 = new Dfloat[N_WE];
  BE0 = new Dfloat[N_WE];
  BW1 = new Dfloat[N_WE];
  BE1 = new Dfloat[N_WE];
  

  BUp0 = new Dfloat[N_UpDown]; 
  BDown0 = new Dfloat[N_UpDown];
  BUp1 = new Dfloat[N_UpDown]; 
  BDown1 = new Dfloat[N_UpDown];

}

MPI_DATA::~MPI_DATA(){
  delete [] BS0;
  delete [] BN0;
  delete [] BW0;
  delete [] BE0;
  delete [] BUp0;
  delete [] BDown0;
  delete [] BS1;
  delete [] BN1;
  delete [] BW1;
  delete [] BE1;
  delete [] BUp1;
  delete [] BDown1;

  sdm = NULL;
}


void MPI_DATA::TRANSFER(char *VarName){

  int indxS,indxD;
  VecI Subindx = sdm->SubIdx();
  VecI SubNum = sdm->SubNum();

  if (sdm->BNorth()){
     
     indxD = Subindx.x + (Subindx.y + 1) * SubNum.x +	\
     Subindx.z * SubNum.x * SubNum.y;
     
     
     //if (sdm->BNorth()) indxD = MPI::PROC_NULL;

     //printf("%d\t%d\n",indxS,indxD);

      sdm->ExpBoundary(BN0,"North",VarName);

      MPI::COMM_WORLD.Sendrecv(BN0,N_SN,MY_MPI_Dfloat,indxD,0,BN1,	\
       		       N_SN,MY_MPI_Dfloat,indxD,0);

      sdm->ImpBoundary(BN1,"North",VarName);

  }
      

  if (sdm->BSouth()){ 

      indxD = Subindx.x + (Subindx.y - 1) * SubNum.x +	\
	Subindx.z * SubNum.x * SubNum.y;
     
      // if (sdm->BSouth() == 0) indxD = MPI::PROC_NULL;

     //printf("%d\t%d\n",indxS,indxD);

      sdm->ExpBoundary(BS0,"South",VarName);

      MPI::COMM_WORLD.Sendrecv(BS0,N_SN,MY_MPI_Dfloat,indxD,0,BS1,	\
			       N_SN,MY_MPI_Dfloat,indxD,0);

      sdm->ImpBoundary(BS1,"South",VarName);
  }
     

    

  if (sdm->BEast()){
    
      indxD = (Subindx.x + 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BEast() == 0) indxD = MPI::PROC_NULL;
     
      sdm->ExpBoundary(BE0,"East",VarName);

      MPI::COMM_WORLD.Sendrecv(BE0,N_WE,MY_MPI_Dfloat,indxD,0,BE1,	\
			       N_WE,MY_MPI_Dfloat,indxD,0);


      sdm->ImpBoundary(BE1,"East",VarName);

  }

   


  if (sdm->BWest()){

      indxD = (Subindx.x - 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BWest() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BW0,"West",VarName);

      MPI::COMM_WORLD.Sendrecv(BW0,N_WE,MY_MPI_Dfloat,indxD,0,BW1,	\
			       N_WE,MY_MPI_Dfloat,indxD,0);

      sdm->ImpBoundary(BW1,"West",VarName);


  }
   


  if (sdm->BUp()){
    
      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z + 1) * SubNum.x * SubNum.y;

      //if (sdm->BUp() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BUp0,"UP",VarName);

      MPI::COMM_WORLD.Sendrecv(BUp0,N_UpDown,MY_MPI_Dfloat,indxD,0,BUp1, \
			       N_UpDown,MY_MPI_Dfloat,indxD,0);

      sdm->ImpBoundary(BUp1,"UP",VarName);

  }
      
    


  if (sdm->BDown()){

      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z -1) * SubNum.x * SubNum.y;

      //if (sdm->BDown() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BDown0,"DOWN",VarName);

      MPI::COMM_WORLD.Sendrecv(BDown0,N_UpDown,MY_MPI_Dfloat,indxD,0,BDown1, \
			       N_UpDown,MY_MPI_Dfloat,indxD,0);

      sdm->ImpBoundary(BDown1,"DOWN",VarName);


  } 

}

void MPI_DATA::Merge(Dfloat *LocVar,int NX,int NY, int NZ,VecI *subi,int rank){

  if (rank == 0) {

    Dfloat *buff,*GlobVar;
    VecI indx;
    int GlobIdx;
    int LocIdx;
    int Nsub = sdm->NumSubDom.x * sdm->NumSubDom.y * sdm->NumSubDom.z; 

    buff = new Dfloat[sdm->SNodeX() * sdm->SNodeY() * sdm->SNodeZ()];
    GlobVar = new Dfloat[NX * NY * NZ];


    for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){

	  indx = sdm->Loc2Glo({i,j,k});
	  GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY; 
	    
	    GlobVar[GlobIdx] = LocVar[sdm->IJK(i,j,k)];


	}
      }
    }


    for (int ii = 1; ii<Nsub;++ii) {


      MPI::COMM_WORLD.Recv(buff,sdm->SNodeT(),MY_MPI_Dfloat,ii,0);
      

      for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
	for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
	  for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){


	    indx.x =  i - KHALO + subi[ii].x * sdm->NodLoc.x;
	    indx.y =  j - KHALO + subi[ii].y * sdm->NodLoc.y;
	    indx.z =  k - KHALO + subi[ii].z * sdm->NodLoc.z;
	    
	    GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY; 
	    GlobVar[GlobIdx] = buff[sdm->IJK(i,j,k)];

	  }
	}
      }


    } // Loop Subdomains


    


  } if (rank >0) {
    

    MPI::COMM_WORLD.Send(LocVar,sdm->SNodeT(),MY_MPI_Dfloat,0,0);


  }


  

}




void MPI_DATA::MergePrint(Dfloat *LocVar,int NX,int NY, int NZ,VecI *subi,int rank, char *name){

  if (rank == 0) {

    Dfloat *buff,*GlobVar;
    VecI indx;
    int GlobIdx;
    int LocIdx;
    int Nsub = sdm->NumSubDom.x * sdm->NumSubDom.y * sdm->NumSubDom.z; 

    buff = new Dfloat[sdm->SNodeX() * sdm->SNodeY() * sdm->SNodeZ()];
    GlobVar = new Dfloat[NX * NY * NZ];


    for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
	for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){

	  indx = sdm->Loc2Glo({i,j,k});
	  GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY; 
	    
	    GlobVar[GlobIdx] = LocVar[sdm->IJK(i,j,k)];


	}
      }
    }


    for (int ii = 1; ii<Nsub;++ii) {


      MPI::COMM_WORLD.Recv(buff,sdm->SNodeT(),MY_MPI_Dfloat,ii,0);
      
      // printf("\n proc %d ; %d \n",Nsub,ii);

      for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
	for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
	  for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){


	    indx.x =  i - KHALO + subi[ii].x * sdm->NodLoc.x;
	    indx.y =  j - KHALO + subi[ii].y * sdm->NodLoc.y;
	    indx.z =  k - KHALO + subi[ii].z * sdm->NodLoc.z;
	    
	    GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY; 
	    GlobVar[GlobIdx] = buff[sdm->IJK(i,j,k)];

	  }
	}
      }


    } // Loop Subdomains


    FILE *R;

    R=fopen(name,"wb");

    for (int k=PML.z;k<NZ;k++){
	for (int j=PML.y;j<NY-PML.y;j++){
	  for (int i=PML.x;i<NX-PML.x;i++){


	     GlobIdx = i + j * NX + k * NX * NY;
	     
	    fwrite(&GlobVar[GlobIdx],sizeof(Dfloat),1,R);

	  }
	}
    }

    fclose(R);

    


    


  } if (rank >0) {
    

    MPI::COMM_WORLD.Send(LocVar,sdm->SNodeT(),MY_MPI_Dfloat,0,0);


  }


  

}

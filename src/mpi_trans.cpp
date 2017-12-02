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

      MPI::COMM_WORLD.Sendrecv(BN0,N_SN,MPI_DOUBLE,indxD,0,BN1,	\
       		       N_SN,MPI_DOUBLE,indxD,0);

      sdm->ImpBoundary(BN1,"North",VarName);

  }
      

  if (sdm->BSouth()){ 

      indxD = Subindx.x + (Subindx.y - 1) * SubNum.x +	\
	Subindx.z * SubNum.x * SubNum.y;
     
      // if (sdm->BSouth() == 0) indxD = MPI::PROC_NULL;

     //printf("%d\t%d\n",indxS,indxD);

      sdm->ExpBoundary(BS0,"South",VarName);

      MPI::COMM_WORLD.Sendrecv(BS0,N_SN,MPI_DOUBLE,indxD,0,BS1,	\
			       N_SN,MPI_DOUBLE,indxD,0);

      sdm->ImpBoundary(BS1,"South",VarName);
  }
     

    

  if (sdm->BEast()){
    
      indxD = (Subindx.x + 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BEast() == 0) indxD = MPI::PROC_NULL;
     
      sdm->ExpBoundary(BE0,"East",VarName);

      MPI::COMM_WORLD.Sendrecv(BE0,N_WE,MPI_DOUBLE,indxD,0,BE1,	\
			       N_WE,MPI_DOUBLE,indxD,0);


      sdm->ImpBoundary(BE1,"East",VarName);

  }

   


  if (sdm->BWest()){

      indxD = (Subindx.x - 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BWest() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BW0,"West",VarName);

      MPI::COMM_WORLD.Sendrecv(BW0,N_WE,MPI_DOUBLE,indxD,0,BW1,	\
			       N_WE,MPI_DOUBLE,indxD,0);

      sdm->ImpBoundary(BW1,"West",VarName);


  }
   


  if (sdm->BUp()){
    
      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z + 1) * SubNum.x * SubNum.y;

      //if (sdm->BUp() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BUp0,"UP",VarName);

      MPI::COMM_WORLD.Sendrecv(BUp0,N_UpDown,MPI_DOUBLE,indxD,0,BUp1, \
			       N_UpDown,MPI_DOUBLE,indxD,0);

      sdm->ImpBoundary(BUp1,"UP",VarName);

  }
      
    


  if (sdm->BDown()){

      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z -1) * SubNum.x * SubNum.y;

      //if (sdm->BDown() == 0) indxD = MPI::PROC_NULL;

      sdm->ExpBoundary(BDown0,"DOWN",VarName);

      MPI::COMM_WORLD.Sendrecv(BDown0,N_UpDown,MPI_DOUBLE,indxD,0,BDown1, \
			       N_UpDown,MPI_DOUBLE,indxD,0);

      sdm->ImpBoundary(BDown1,"DOWN",VarName);


  } 

}

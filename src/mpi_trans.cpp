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
  BS0 = new Dfloat[N_SN * 6]; 
  BN0 = new Dfloat[N_SN * 6];
  BS1 = new Dfloat[N_SN * 6]; 
  BN1 = new Dfloat[N_SN * 6];
  
  BW0 = new Dfloat[N_WE * 6];
  BE0 = new Dfloat[N_WE * 6];
  BW1 = new Dfloat[N_WE * 6];
  BE1 = new Dfloat[N_WE * 6];

  BUp0 = new Dfloat[N_UpDown * 6]; 
  BDown0 = new Dfloat[N_UpDown * 6];
  BUp1 = new Dfloat[N_UpDown * 6]; 
  BDown1 = new Dfloat[N_UpDown * 6];

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


void MPI_DATA::TRANSFER(int Var){

  int indxS,indxD;
  VecI Subindx = sdm->SubIdx();
  VecI SubNum = sdm->SubNum();
  int nvar;
  int nchar;

  	if ( Var == 1 ) {
	// VELOCITIES
	nchar = 0;
	nvar = 6;

	} else if (Var == 2){
	// STRESESS
	nchar = 6;
	nvar = 6;

	}
  

  if (sdm->BNorth()){
     
     indxD = Subindx.x + (Subindx.y + 1) * SubNum.x +	\
     Subindx.z * SubNum.x * SubNum.y;
     
     
     //if (sdm->BNorth()) indxD = MPI::PROC_NULL;

     //printf("%d\t%d\n",indxS,indxD);


	for (int ivar = 0 ;ivar<nvar;++ivar){
      sdm->ExpBoundary((BN0 + ivar * N_SN),"North",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BN0,N_SN * nvar,MY_MPI_Dfloat,indxD,0,BN1,	\
       		       N_SN * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BN1 + ivar * N_SN),"North",VarName[ivar + nchar]);
	}

  }
      

  if (sdm->BSouth()){ 

      indxD = Subindx.x + (Subindx.y - 1) * SubNum.x +	\
	Subindx.z * SubNum.x * SubNum.y;
     
      // if (sdm->BSouth() == 0) indxD = MPI::PROC_NULL;

     //printf("%d\t%d\n",indxS,indxD);


	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ExpBoundary((BS0 + ivar * N_SN),"South",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BS0,N_SN * nvar,MY_MPI_Dfloat,indxD,0,BS1,	\
			       N_SN * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BS1 + ivar * N_SN),"South",VarName[ivar + nchar]);
	}

  }
     

    

  if (sdm->BEast()){
    
      indxD = (Subindx.x + 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BEast() == 0) indxD = MPI::PROC_NULL;


	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ExpBoundary((BE0 + ivar * N_WE),"East",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BE0,N_WE * nvar,MY_MPI_Dfloat,indxD,0,BE1,	\
			       N_WE * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BE1 + ivar * N_WE),"East",VarName[ivar + nchar]);
	}

  }

   


  if (sdm->BWest()){

      indxD = (Subindx.x - 1) + Subindx.y * SubNum.x +	\
       Subindx.z * SubNum.x * SubNum.y;

      //if (sdm->BWest() == 0) indxD = MPI::PROC_NULL;


	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ExpBoundary((BW0 + ivar * N_WE),"West",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BW0,N_WE * nvar,MY_MPI_Dfloat,indxD,0,BW1,	\
			       N_WE * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BW1 + ivar * N_WE),"West",VarName[ivar + nchar]);
	}

  }
   


  if (sdm->BUp()){
    
      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z + 1) * SubNum.x * SubNum.y;

      //if (sdm->BUp() == 0) indxD = MPI::PROC_NULL;



	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ExpBoundary((BUp0 + ivar * N_UpDown),"UP",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BUp0,N_UpDown * nvar,MY_MPI_Dfloat,indxD,0,BUp1, \
			       N_UpDown * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BUp1 + ivar * N_UpDown),"UP",VarName[ivar + nchar]);
	}


  }
      
    


  if (sdm->BDown()){

      indxD = Subindx.x + Subindx.y * SubNum.x +	\
	(Subindx.z -1) * SubNum.x * SubNum.y;

      //if (sdm->BDown() == 0) indxD = MPI::PROC_NULL;


	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ExpBoundary((BDown0 + ivar * N_UpDown),"DOWN",VarName[ivar + nchar]);
	}

      MPI_Sendrecv(BDown0,N_UpDown * nvar,MY_MPI_Dfloat,indxD,0,BDown1, \
			       N_UpDown * nvar,MY_MPI_Dfloat,indxD,0,MPI_COMM_WORLD,&status);

	for (int ivar = 0;ivar<nvar;++ivar){
      sdm->ImpBoundary((BDown1 + ivar * N_UpDown),"DOWN",VarName[ivar + nchar]);
	}

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


      MPI_Recv(buff,sdm->SNodeT(),MY_MPI_Dfloat,ii,0,MPI_COMM_WORLD,&status);
      

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
    

    MPI_Send(LocVar,sdm->SNodeT(),MY_MPI_Dfloat,0,0,MPI_COMM_WORLD);


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


      MPI_Recv(buff,sdm->SNodeT(),MY_MPI_Dfloat,ii,0,MPI_COMM_WORLD,&status);
      
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
    

    MPI_Send(LocVar,sdm->SNodeT(),MY_MPI_Dfloat,0,0,MPI_COMM_WORLD);


  }

}


 void MPI_DATA::KernPrint(Dfloat *LocVar,int NX,int NY, int NZ,VecI *subi,int rank, char *name){

 if (rank == 0) {

    Dfloat *buff,*kr,*kp,*ks,*kpa,*kpb;
    VecI indx;
    int GlobIdx;
    int LocIdx;
    int Nsub = sdm->NumSubDom.x * sdm->NumSubDom.y * sdm->NumSubDom.z;

    buff = new Dfloat[sdm->SNodeT() * 5];
    kr = new Dfloat[(NX * NY * NZ)];
    kp = new Dfloat[(NX * NY * NZ)];
    ks = new Dfloat[(NX * NY * NZ)];
    kpa = new Dfloat[(NX * NY * NZ)];
    kpb = new Dfloat[(NX * NY * NZ)];

    


    for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
      for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
        for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){

          indx = sdm->Loc2Glo({i,j,k});
          GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY;

            kr[GlobIdx] = LocVar[sdm->IJK(i,j,k)];
	    kp[GlobIdx] = LocVar[sdm->IJK(i,j,k) + sdm->SNodeT()];
	    ks[GlobIdx] = LocVar[sdm->IJK(i,j,k) + 2 * sdm->SNodeT()];
	    kpa[GlobIdx] = LocVar[sdm->IJK(i,j,k) + 3 * sdm->SNodeT()];
	    kpb[GlobIdx] = LocVar[sdm->IJK(i,j,k) + 4 * sdm->SNodeT()];


        }
      }
    }




 for (int ii = 1; ii<Nsub;++ii) {


      MPI_Recv(buff,sdm->SNodeT() * 5,MY_MPI_Dfloat,ii,0,MPI_COMM_WORLD,&status);

      // printf("\n proc %d ; %d \n",Nsub,ii);

      for (int k=KHALO;k<sdm->SNodeZ()-KHALO;k++){
        for (int j=KHALO;j<sdm->SNodeY()-KHALO;j++){
          for (int i=KHALO;i<sdm->SNodeX()-KHALO;i++){


            indx.x =  i - KHALO + subi[ii].x * sdm->NodLoc.x;
            indx.y =  j - KHALO + subi[ii].y * sdm->NodLoc.y;
            indx.z =  k - KHALO + subi[ii].z * sdm->NodLoc.z;

            GlobIdx = indx.x + indx.y * NX + indx.z * NX * NY;
            
	    kr[GlobIdx] = buff[sdm->IJK(i,j,k)];
            kp[GlobIdx] = buff[sdm->IJK(i,j,k) + sdm->SNodeT()];
            ks[GlobIdx] = buff[sdm->IJK(i,j,k) + 2 * sdm->SNodeT()];
            kpa[GlobIdx] = buff[sdm->IJK(i,j,k) + 3 * sdm->SNodeT()];
            kpb[GlobIdx] = buff[sdm->IJK(i,j,k) + 4 * sdm->SNodeT()];

          }
        }
      }


  }




 FILE *R,*S,*P,*PA,*PB;
 char rn[100],pn[100],sn[100],pan[100],pbn[100];
 Dfloat *kr_buff,*kp_buff,*ks_buff,*kpa_buff,*kpb_buff;
 int size_buff = (NZ - PML.z) * (NX - 2 * PML.x) * (NY - 2 * PML.y);
 int idx;

 kr_buff = new Dfloat[size_buff];
 kp_buff = new Dfloat[size_buff];
 ks_buff = new Dfloat[size_buff];
 kpa_buff = new Dfloat[size_buff];
 kpb_buff = new Dfloat[size_buff];

 snprintf(rn,90,"%sKRHO.bin",name);
 snprintf(pn,90,"%sKVP.bin",name);
 snprintf(sn,90,"%sKVS.bin",name);
 snprintf(pan,90,"%sPcondA.bin",name);
 snprintf(pbn,90,"%sPcondB.bin",name);

    R=fopen(rn,"wb");
    P=fopen(pn,"wb");
    S=fopen(sn,"wb");
    PA=fopen(pan,"wb");
    PB=fopen(pbn,"wb");

    for (int k=PML.z;k<NZ;k++){
        for (int j=PML.y;j<NY-PML.y;j++){
          for (int i=PML.x;i<NX-PML.x;i++){


             GlobIdx = i + j * NX + k * NX * NY;
    	     idx =  (i - PML.x) + (j - PML.y) * (NX - 2 * PML.x) + \
	       (k - PML.z) * (NX - 2 * PML.x) * (NY - 2 * PML.y);
	     
             kr_buff[idx] = kr[GlobIdx];
    	     kp_buff[idx] = kp[GlobIdx];
    	     ks_buff[idx] = ks[GlobIdx];
    	     kpa_buff[idx] = kpa[GlobIdx];
    	     kpb_buff[idx] = kpb[GlobIdx];

	    //  fwrite(&kr[GlobIdx],sizeof(Dfloat),1,R);
	    // fwrite(&kp[GlobIdx],sizeof(Dfloat),1,P);
	    // fwrite(&ks[GlobIdx],sizeof(Dfloat),1,S);
	    // fwrite(&kpa[GlobIdx],sizeof(Dfloat),1,PA);
	    // fwrite(&kpb[GlobIdx],sizeof(Dfloat),1,PB);

          }
        }
    }

     fwrite(kr_buff,sizeof(Dfloat),size_buff,R);
     fwrite(kp_buff,sizeof(Dfloat),size_buff,P);
     fwrite(ks_buff,sizeof(Dfloat),size_buff,S);
     fwrite(kpa_buff,sizeof(Dfloat),size_buff,PA);
     fwrite(kpb_buff,sizeof(Dfloat),size_buff,PB);
  

    fclose(R);
    fclose(P);
    fclose(S);
    fclose(PA);
    fclose(PB);

    delete [] kr;
    delete [] kp;
    delete [] ks;
    delete [] kpa;
    delete [] kpb;
    delete [] buff;

    delete [] kr_buff;
    delete [] kp_buff;
    delete [] ks_buff;
    delete [] kpa_buff;
    delete [] kpb_buff;

}

if (rank >0) {


    MPI_Send(LocVar,sdm->SNodeT() * 5,MY_MPI_Dfloat,0,0,MPI_COMM_WORLD);


  }


}

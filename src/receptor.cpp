#include "receptor.hpp"


receptor::receptor(geometry3D *domain, std::string nFile,int nrecep,int in_nt) {

  nt = in_nt;
  GDomain  = domain;
  FileS = nFile;
  nr = nrecep;
  xcoord = new Dfloat[nr];
  ycoord = new Dfloat[nr];
  zcoord = new Dfloat[nr];
  pos_recep = new VecI[nr];

  vx_ad = new Dfloat[nt * nr];
  vy_ad = new Dfloat[nt * nr];
  vz_ad = new Dfloat[nt * nr];
  
  nameStation = new std::string[nr];
  RX = new std::fstream[nr];
  RY = new std::fstream[nr];
  RZ = new std::fstream[nr];
  
  //std::cout<<FileS.c_str()<<std::endl;
  
  R.open(FileS.c_str());
  int npos1 = 0;
  int npos2 = 0;
  std::getline(R,line);
  
  // std::cout<<line.c_str()<<std::endl;
  for (int i=0; i<nr; ++i){
    //fgets(data,200,R);

    std::getline(R,line);
    //std::cout<<line<<std::endl;
    npos2 = line.find(",",npos1);
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    nameStation[i] = line.substr(npos1,npos2-npos1);
    npos1 = npos2 + 1;
    //std::cout<<" X:"<<xcoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    

    npos2 = line.find(",",npos1);
    // std::cout<<npos1<<" :"<<npos2<<std::endl;
    xcoord[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" Y:"<<ycoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
   
    npos2 = line.find(",",npos1);
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    ycoord[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" Z:"<<zcoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    
    npos2 = line.find(",",npos1);
    zcoord[i] = std::stof(line.substr(npos1));
    //std::cout<<" azimuth:"<<azimuth[i]<<std::endl;

   npos1 = 0;
   npos2 = 0;

    pos_recep[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]});
  
   //std::cout<<pos_recep[i].x<<pos_recep[i].y<<pos_recep[i].z<<"Receptor"<<std::endl;  
  }
  

  
    
}

void receptor::FileOpen(int i, int PROPAGATION){

  if (!PROPAGATION){
  RX[i].open("DATA/"+nameStation[i] + "-VX.bin",std::ios::binary | std::ios::trunc | std::ios::out);
  RY[i].open("DATA/"+nameStation[i] + "-VY.bin",std::ios::binary | std::ios::trunc | std::ios::out);
  RZ[i].open("DATA/"+nameStation[i] + "-VZ.bin",std::ios::binary | std::ios::trunc | std::ios::out);
  }

  if (PROPAGATION){
    
  RX[i].open("DATA/"+nameStation[i] + "-ADJX.bin",std::ios::binary | std::ios::in);
  RY[i].open("DATA/"+nameStation[i] + "-ADJY.bin",std::ios::binary | std::ios::in);
  RZ[i].open("DATA/"+nameStation[i] + "-ADJZ.bin",std::ios::binary | std::ios::in);

  }
  

}

void receptor::FileClose(int i){

  RX[i].close();
  RY[i].close();
  RZ[i].close();
  
  

}

void receptor::WriteFile(int i, Dfloat vx, Dfloat vy,Dfloat vz){

  RX[i].write( (char*)&vx, sizeof(Dfloat));
  RY[i].write( (char*)&vy, sizeof(Dfloat));
  RZ[i].write( (char*)&vz, sizeof(Dfloat));

}


void receptor::LoadFile(int i){


  for (int it = 0; it<nt;++it){
    RX[i].read( (char*)&vx_ad[it + nt * i], sizeof(Dfloat));
    RY[i].read( (char*)&vy_ad[it + nt * i], sizeof(Dfloat));
    RZ[i].read( (char*)&vz_ad[it + nt * i], sizeof(Dfloat));
  }


  RX[i].close();
  RY[i].close();
  RZ[i].close();

  
}





receptor::~receptor(){

  delete [] RX;
  delete [] RY;
  delete [] RZ;
  delete [] nameStation;
  delete [] xcoord;
  delete [] ycoord;
  delete [] zcoord;
  delete [] vx_ad;
  delete [] vy_ad;
  delete [] vz_ad;
  

  GDomain = NULL;

  R.close();
  
}




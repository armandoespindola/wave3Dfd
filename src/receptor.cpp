#include "receptor.hpp"


receptor::receptor(geometry3D *domain, std::string nFile,int nrecep) {

  GDomain  = domain;
  FileS = nFile;
  nr = nrecep;
  xcoord = new Dfloat[nr];
  ycoord = new Dfloat[nr];
  zcoord = new Dfloat[nr];
  pos_recep = new VecI[nr];
  
  nameStation = new std::string[nr];
  RX = new std::ofstream[nr];
  RY = new std::ofstream[nr];
  RZ = new std::ofstream[nr];
  
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

    pos_recep[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]});
    
  }
  

  
    
}

void receptor::FileOpen(int i){

  RX[i].open("../data/"+nameStation[i] + "-VX.bin",std::ios::binary | std::ios::trunc);
  RY[i].open("../data/"+nameStation[i] + "-VY.bin",std::ios::binary | std::ios::trunc);
  RZ[i].open("../data/"+nameStation[i] + "-VZ.bin",std::ios::binary | std::ios::trunc);
  
  

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





receptor::~receptor(){

  delete [] RX;
  delete [] RY;
  delete [] RZ;
  delete [] nameStation;
  delete [] xcoord;
  delete [] ycoord;
  delete [] zcoord;

  R.close();
  
}




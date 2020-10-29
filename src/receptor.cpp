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
  pos_vx = new VecI[nr];
  pos_vy = new VecI[nr];
  pos_vz = new VecI[nr];
  pos_sii = new VecI[nr];
  pos_sxy = new VecI[nr];
  pos_sxz = new VecI[nr];
  pos_syz = new VecI[nr];

  vx_ad = new Dfloat[nt * nr];
  vy_ad = new Dfloat[nt * nr];
  vz_ad = new Dfloat[nt * nr];


  for (int i=0; i<nr*nt; ++i){
    vx_ad[i] = 0.0;
    vy_ad[i] = 0.0;
    vz_ad[i] = 0.0;
  }
  
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

   pos_recep[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,1,1});
   pos_vx[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,0,0});
   pos_vy[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,1,0});
   pos_vz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,0,1});
   pos_sii[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,0,0});
   pos_sxy[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,1,0});
   pos_sxz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,0,1});
   pos_syz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,1,1});
  
   // std::cout<<i<<" "<<pos_vx[i].x<<" "<<pos_vx[i].y<<" "<<pos_vx[i].z<<"Receptor"<<std::endl;
   // std::cout<<i<<" "<<pos_vy[i].x<<" "<<pos_vy[i].y<<" "<<pos_vy[i].z<<"Receptor"<<std::endl;
   // std::cout<<i<<" "<<pos_vz[i].x<<" "<<pos_vz[i].y<<" "<<pos_vz[i].z<<"Receptor"<<std::endl; 
  }
  

  
    
}

void receptor::FileOpen(int i, int PROPAGATION,char *NameVar){

  if (!PROPAGATION){

    if  (strcmp("VX",NameVar) == 0){
  RX[i].open("DATA/"+nameStation[i] + "-VX.bin",std::ios::binary | std::ios::trunc | std::ios::out);
    } else  if (strcmp("VY",NameVar) == 0){
  RY[i].open("DATA/"+nameStation[i] + "-VY.bin",std::ios::binary | std::ios::trunc | std::ios::out);
    } else if  (strcmp("VZ",NameVar) == 0){
  RZ[i].open("DATA/"+nameStation[i] + "-VZ.bin",std::ios::binary | std::ios::trunc | std::ios::out);
    }

  }

  if (PROPAGATION){
    
  // RX[i].open("DATA/"+nameStation[i] + "-VX.bin",std::ios::binary | std::ios::in);
  // RY[i].open("DATA/"+nameStation[i] + "-VY.bin",std::ios::binary | std::ios::in);
  // RZ[i].open("DATA/"+nameStation[i] + "-VZ.bin",std::ios::binary | std::ios::in);

    if  (strcmp("VX",NameVar) == 0){
  RX[i].open("DATA/"+nameStation[i] + "-ADJX.bin",std::ios::binary | std::ios::in);
    }else if  (strcmp("VY",NameVar) == 0){
  RY[i].open("DATA/"+nameStation[i] + "-ADJY.bin",std::ios::binary | std::ios::in);
    } else if  (strcmp("VZ",NameVar) == 0){
  RZ[i].open("DATA/"+nameStation[i] + "-ADJZ.bin",std::ios::binary | std::ios::in);
    }

  }
  

  }

  void receptor::FileClose(int i,char *NameVar){

  
     if  (strcmp("VX",NameVar) == 0){
      RX[i].close();
    }else if  (strcmp("VY",NameVar) == 0){
       RY[i].close();
    }else if  (strcmp("VZ",NameVar) == 0){
       RZ[i].close();
    }
  

  
  

}

  void receptor::WriteFile(int i, Dfloat *var,char *NameVar){


     if  (strcmp("VX",NameVar) == 0){
       RX[i].write( (char*)&var[0], sizeof(Dfloat) * nt);
    }else if  (strcmp("VY",NameVar) == 0){
      RY[i].write( (char*)&var[0], sizeof(Dfloat) * nt );
    }else if  (strcmp("VZ",NameVar) == 0){
     RZ[i].write( (char*)&var[0], sizeof(Dfloat) * nt); 
    }

}


  void receptor::LoadFile(int i,char *NameVar){

    if  (strcmp("VX",NameVar) == 0){
       RX[i].read( (char*)&vx_ad[nt * i], sizeof(Dfloat) * nt);
       RX[i].close();
      
    }else if  (strcmp("VY",NameVar) == 0){
       RY[i].read( (char*)&vy_ad[nt * i], sizeof(Dfloat) * nt);
        RY[i].close();
    }else if  (strcmp("VZ",NameVar) == 0){
      RZ[i].read( (char*)&vz_ad[nt * i], sizeof(Dfloat) * nt);
       RZ[i].close();
    }
    
  //for (int it = 0; it<nt;++it){
   
  
   
    //  }


  
 

  
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
  delete [] pos_recep;
  delete [] pos_vx;
  delete [] pos_vy;
  delete [] pos_vz;
  delete [] pos_sii;
  delete [] pos_sxy;
  delete [] pos_sxz;
  delete [] pos_syz;
  

  GDomain = NULL;

  R.close();
  
}


void receptor::PrintInf(){

  Dfloat xm,ym,zm;
      
  std::cout<<""<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<"######################"<<std::endl;
  std::cout<<"RECEPTORS INFORMATION"<<std::endl;
  std::cout<<"######################"<<std::endl;
  
  
  for (int i=0; i<nr; ++i){

    std::cout<<""<<std::endl;
    std::cout<<"##########################"<<std::endl;
    std::cout<<" STATION # "<<i + 1<<" : "<<nameStation[i]<<std::endl;
    
    std::cout<<" X: "<<xcoord[i]<<std::endl;
    std::cout<<" Y: "<<ycoord[i]<<std::endl;
    std::cout<<" Z: "<<zcoord[i]<<std::endl;
    std::cout<<" Node_X(Mesh Location) : "<<pos_recep[i].x<<std::endl;
    std::cout<<" Node_Y(Mesh Location) : "<<pos_recep[i].y<<std::endl;
    std::cout<<" Node_Z(Mesh Location) : "<<GDomain->HALO_NodeZ()-pos_recep[i].z-1<<std::endl;

    xm = GDomain->CoorXHalf(pos_recep[i].x);
    ym = GDomain->CoorYHalf(pos_recep[i].y);
    zm = GDomain->CoorZHalf(GDomain->HALO_NodeZ() - 1 -pos_recep[i].z);

    std::cout<<" X(Mesh Location) : "<<xm<<std::endl;
    std::cout<<" Y(Mesh Location) : "<<ym<<std::endl;
    std::cout<<" Z(Mesh Location) : "<<zm<<std::endl;

    std::cout<<" X(Error Location) : "<<abs(xm - xcoord[i]) <<std::endl;
    std::cout<<" Y(Error Location) : "<<abs(ym - ycoord[i])<<std::endl;
    std::cout<<" Z(Error Location) : "<<abs(zm - zcoord[i])<<std::endl;
    
  }
  
  
}




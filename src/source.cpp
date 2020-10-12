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

#include "source.hpp"


source::source(geometry3D *domain, std::string nFile,int nsource) {

  GDomain  = domain;
  FileS = nFile;
  ns = nsource;
  Mxx = new Dfloat[ns];
  Myy = new Dfloat[ns];
  Mzz = new Dfloat[ns];
  Mxy = new Dfloat[ns];
  Mxz = new Dfloat[ns];
  Myz = new Dfloat[ns];
  nshift = new int[ns];
  xcoord = new Dfloat[ns];
  ycoord = new Dfloat[ns];
  zcoord = new Dfloat[ns];

  pos_vx = new VecI[ns];
  pos_vy = new VecI[ns];
  pos_vz = new VecI[ns];
  pos_sii = new VecI[ns];
  pos_sxy = new VecI[ns];
  pos_sxz = new VecI[ns];
  pos_syz = new VecI[ns];

  pos_src = new VecI[ns];
  src_r_half = new Dfloat[ns];
  src_r = new Dfloat[ns];
  src_flag_r = new int[ns];
    
  coor_vx = new VecF[ns];
  coor_vy = new VecF[ns];
  coor_vz = new VecF[ns];
  coor_sii = new VecF[ns];
  coor_sxy = new VecF[ns];
  coor_sxz = new VecF[ns];
  coor_syz = new VecF[ns];
  
  M0 = new Dfloat[ns];
  d_t0 = new Dfloat[ns];
  strike = new Dfloat[ns];
  dip = new Dfloat[ns];
  slip = new Dfloat[ns];
  azimuth = new Dfloat[ns];
  nameSource = new std::string[ns];

  //std::cout<<FileS.c_str()<<std::endl;
  
  R.open(FileS.c_str());
  int npos1 = 0;
  int npos2 = 0;
  std::getline(R,line);
  
  // std::cout<<line.c_str()<<std::endl;
  for (int i=0; i<ns; ++i){
    //fgets(data,200,R);

    std::getline(R,line);
    //std::cout<<line<<std::endl;
    npos2 = line.find(",",npos1);
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    nameSource[i] = line.substr(npos1,npos2-npos1);
    npos1 = npos2 + 1;

    npos2 = line.find(",",npos1);
    xcoord[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" X:"<<xcoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    

    npos2 = line.find(",",npos1);
    // std::cout<<npos1<<" :"<<npos2<<std::endl;
    ycoord[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" Y:"<<ycoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
   
    npos2 = line.find(",",npos1);
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    zcoord[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" Z:"<<zcoord[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    
    npos2 = line.find(",",npos1);
    //std::cout<<npos1<<" :"<<npos2<<std::endl;
    M0[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" M0:"<<M0[i];
    //std::cout<<npos1<<" :"<<npos2<<std::endl;

    
    
    npos2 = line.find(",",npos1);
    d_t0[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" d_t0:"<<d_t0[i];

    npos2 = line.find(",",npos1);
    Mxx[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" strike:"<<strike[i];

    npos2 = line.find(",",npos1);
    Mxy[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" strike:"<<strike[i];

    npos2 = line.find(",",npos1);
    Mxz[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" strike:"<<strike[i];

    npos2 = line.find(",",npos1);
    Myy[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" dip:"<<dip[i];

    npos2 = line.find(",",npos1);
    Myz[i] = std::stof(line.substr(npos1,npos2-npos1));
    npos1 = npos2 + 1;
    //std::cout<<" slip:"<<slip[i];
    
    npos2 = line.find(",",npos1);
    Mzz[i] = std::stof(line.substr(npos1));
    //std::cout<<" azimuth:"<<azimuth[i]<<std::endl;

    npos1 = 0;
    npos2 = 0;


    pos_vx[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,0,0});
    coor_vx[i] = {GDomain->CoorX(pos_vx[i].x),GDomain->CoorY(pos_vx[i].y),\
                  GDomain->CoorZ(GDomain->HALO_NodeZ() - 1 - pos_vx[i].z)};


    pos_vy[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,1,0});
    coor_vy[i] = {GDomain->CoorXHalf(pos_vy[i].x),GDomain->CoorYHalf(pos_vy[i].y),\
                  GDomain->CoorZ(GDomain->HALO_NodeZ() - 1 - pos_vy[i].z)};

    pos_vz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,0,1});
    coor_vz[i] = {GDomain->CoorXHalf(pos_vz[i].x),GDomain->CoorY(pos_vz[i].y),\
                  GDomain->CoorZHalf(GDomain->HALO_NodeZ() - 1 - pos_vz[i].z)};

    pos_sii[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,0,0});
    coor_sii[i] = {GDomain->CoorXHalf(pos_sii[i].x),GDomain->CoorY(pos_sii[i].y),\
                  GDomain->CoorZ(GDomain->HALO_NodeZ() - 1 - pos_sii[i].z)};

    pos_sxy[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,1,0});
    coor_sxy[i] = {GDomain->CoorX(pos_sxy[i].x),GDomain->CoorYHalf(pos_sxy[i].y),       \
                  GDomain->CoorZ(GDomain->HALO_NodeZ() - 1 - pos_sxy[i].z)};

    pos_sxz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{0,0,1});
    coor_sxz[i] = {GDomain->CoorX(pos_sxz[i].x),GDomain->CoorY(pos_sxz[i].y),\
                  GDomain->CoorZHalf(GDomain->HALO_NodeZ() - 1 - pos_sxz[i].z)};

    pos_syz[i] = GDomain->FindNode({xcoord[i],ycoord[i],zcoord[i]},{1,1,1});
    coor_syz[i] = {GDomain->CoorXHalf(pos_syz[i].x),GDomain->CoorYHalf(pos_syz[i].y),\
                  GDomain->CoorZHalf(GDomain->HALO_NodeZ() - 1 - pos_syz[i].z)};


    src_r_half[i] = sqrt(pow(coor_syz[i].x - xcoord[i],2.0) + \
			 pow(coor_syz[i].y - ycoord[i],2.0) + \
			 pow(coor_syz[i].z - zcoord[i],2.0));
    
    src_r[i] = sqrt(pow(coor_vx[i].x - xcoord[i],2.0) + \
		    pow(coor_vx[i].y - ycoord[i],2.0) +	\
		    pow(coor_vx[i].z - zcoord[i],2.0));


    // Selecting best candidate for source position

    //std::cout<<src_r_half[i]<<" "<<src_r[i]<<std::endl;
    if (src_r_half[i] <= src_r[i]){
      src_flag_r[i] = -1;
      pos_src[i] = pos_syz[i];
      //  std::cout<<" f "<<src_flag_r[i]<<std::endl;
    } else {
      src_flag_r[i] = 1;
      pos_src[i] = pos_vx[i];
      // std::cout<<" f "<<src_flag_r[i]<<std::endl;
	}
      
     
    nshift[i] = GDomain->HALO_NodeZ()-pos_src[i].z-1;

    // std::cout<<"sii "<<coor_sii[i].x<<" "<<coor_sii[i].y<<" "<<coor_sii[i].z<<std::endl;
    // std::cout<<"sxy "<<coor_sxy[i].x<<" "<<coor_sxy[i].y<<" "<<coor_sxy[i].z<<std::endl;
    // std::cout<<"sxz "<<coor_sxz[i].x<<" "<<coor_sxz[i].y<<" "<<coor_sxz[i].z<<std::endl;
    // std::cout<<"syz "<<coor_syz[i].x<<" "<<coor_syz[i].y<<" "<<coor_syz[i].z<<std::endl;
    // std::cout<<"pos_syz "<<pos_syz[i].x<<" "<<pos_syz[i].y<<" "<<pos_syz[i].z<<std::endl;

//   std::cout<<pos_src[i].x<<pos_src[i].y<<pos_src[i].z<<"SOURCE"<<std::endl;
   //  std::cout<<pos_vx[i].x<<" "<<pos_vx[i].y<<" "<<pos_vx[i].z<<"Source"<<std::endl;
   //  std::cout<<pos_vy[i].x<<" "<<pos_vy[i].y<<" "<<pos_vy[i].z<<"Source"<<std::endl;
   // std::cout<<pos_vz[i].x<<" "<<pos_vz[i].y<<" "<<pos_vz[i].z<<"Source"<<std::endl;
    
  }
  
  
  //  smoment();

  
    
}



source::~source(){

  delete [] Mxx;
  delete [] Myy;
  delete [] Mzz;
  delete [] Mxy;
  delete [] Mxz;
  delete [] Myz;
  delete [] M0;
  delete [] d_t0;
  delete [] xcoord;
  delete [] ycoord;
  delete [] zcoord;
  delete [] strike;
  delete [] dip;
  delete [] slip;
  delete [] azimuth;
  delete [] nshift;
  delete [] pos_vx;
  delete [] pos_vy;
  delete [] pos_vz;
  delete [] pos_sii;
  delete [] pos_sxy;
  delete [] pos_sxz;
  delete [] pos_syz;

  delete [] pos_src;
  delete [] src_r_half;
  delete [] src_r;
  delete [] src_flag_r;

  delete [] coor_vx;
  delete [] coor_vy;
  delete [] coor_vz;
  delete [] coor_sii;
  delete [] coor_sxy;
  delete [] coor_sxz;
  delete [] coor_syz;
  
  R.close();
  
}



// Seismic Moment Calculation Aki and Richards (2002) Convention
void source::smoment(){

  Dfloat grd2rad = pi / 180.0 ;
  Dfloat str,d,rake;
 
  for (int i=0;i<ns;i++){
    
    str = grd2rad * (strike[i] - azimuth[i]);
    d = grd2rad * dip[i];
    rake = grd2rad * slip[i];

    
      Mxx[i] =sin(d)*cos(rake)*sin(2.0*str) + sin(2.0*d)*sin(rake)*pow(sin(str),2.0);
      Myy[i] =sin(d)*cos(rake)*sin(2.0*str) - sin(2.0*d)*sin(rake)*pow(cos(str),2.0);
      //Mzz[i] = Mxx[i] + Myy[i];
      Mxy[i] = sin(d)*cos(rake)*cos(2.0*str) + 0.5*sin(2.0*d)*sin(rake)*sin(2.0*str);
      Mxz[i] = cos(d)*cos(rake)*cos(str) + cos(2.0*d)*sin(rake)*sin(str) ;
      Myz[i] = cos(d)*cos(rake)*sin(str) - cos(2.0*d)*sin(rake)*cos(str) ;  
      Mxx[i] *= -M0[i] ; 
      Myy[i] *=  M0[i] ; 
      Mzz[i] = -(Mxx[i] + Myy[i]); 
      Mxy[i] *= M0[i] ; 
      Mxz[i] *= - M0[i] ; 
      Myz[i] *= - M0[i] ; 

//std::cout<<"M11,M22,M33,M12,M13,M23: "<<std::endl;
//std::cout<<Mxx[i] <<Myy[i] <<Mzz[i] <<Mxy[i] <<Mxz[i] <<Myz[i] <<std::endl;
  }

  
}






Dfloat source::sourceType(Dfloat t0, Dfloat f0 , int itime,Dfloat dt, int T_SRC) {

  Dfloat src,a_fu,amp,time;

  time = itime  * dt ; 
  a_fu= pow (pi*f0,2.0);
  src = 0.0;

  // GAUSSIAN 
  
  if (T_SRC==0){ 
    src = exp(-a_fu * pow(time - t0,2.0)) * pow( a_fu / pi , 0.5);
    //src = exp(-a_fu * pow(time - t0,2.0));
  }
  
  // FIRST DERIVATIVE OF A GAUSSIAN
  
  if (T_SRC==1){
    src = exp(-a_fu * pow(time - t0,2.0)) * pow( a_fu / pi , 0.5);
    src *= -2.0 * a_fu * (time - t0); 
  }

  // SECOND DERIVATIVE OF A GAUSSIAN (RICKER PULSE)
  
  if (T_SRC==2){
    src = (1.0 - 2.0 * a_fu * pow((time - t0),2.0)) * exp(-a_fu * pow((time - t0),2.0)) ;
  }

  // HEAVISIDE STEP FUNCTION
  
  if (T_SRC==3){
    if (t0 <= 0.0) {
      if (itime >= 1){
	src = 1.0;
      } else {
      src = 0.0;
      }
    }
  }

  
  return src;

}


void source::PrintInf(){

  Dfloat xm,ym,zm;
      
  std::cout<<""<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<"##################"<<std::endl;
  std::cout<<"SOURCE INFORMATION"<<std::endl;
  std::cout<<"##################"<<std::endl;
  
  
  for (int i=0; i<ns; ++i){

    std::cout<<""<<std::endl;
    std::cout<<"##################"<<std::endl;
    std::cout<<"Source #: "<<i + 1<<" : "<<nameSource[i]<<std::endl;
    
    std::cout<<" X: "<<xcoord[i]<<std::endl;
    std::cout<<" Y: "<<ycoord[i]<<std::endl;
    std::cout<<" Z: "<<zcoord[i]<<std::endl;
    std::cout<<" M0: "<<M0[i]<<std::endl;
    std::cout<<" MXX,MXY,MXZ,MYY,MYZ,MZZ: "<<std::endl;
    std::cout<<" "<<Mxx[i]<<", "<<Mxy[i]<<", "	\
	     <<Mxz[i]<<", "<<Myy[i]<<", "\
	     <<Myz[i]<<", "<<Mzz[i] <<std::endl;

    std::cout<<" Node_X(Mesh Location) : "<<pos_src[i].x<<std::endl;
    std::cout<<" Node_Y(Mesh Location) : "<<pos_src[i].y<<std::endl;
    std::cout<<" Node_Z(Mesh Location) : "<<GDomain->HALO_NodeZ()-pos_src[i].z-1<<std::endl;


    xm = GDomain->CoorXHalf(pos_vx[i].x);
    ym = GDomain->CoorYHalf(pos_vx[i].y);
    zm = GDomain->CoorZHalf(GDomain->HALO_NodeZ() - 1 -pos_vx[i].z);

    std::cout<<" X(Mesh Location) : "<<xm<<std::endl;
    std::cout<<" Y(Mesh Location) : "<<ym<<std::endl;
    std::cout<<" Z(Mesh Location) : "<<zm<<std::endl;

    std::cout<<" X(Error Location) : "<<abs(xm - xcoord[i]) <<std::endl;
    std::cout<<" Y(Error Location) : "<<abs(ym - ycoord[i])<<std::endl;
    std::cout<<" Z(Error Location) : "<<abs(zm - zcoord[i])<<std::endl;
    
  }

  
}


void source::w_sinc(int freeSurf,int is){
  Dfloat arg;
  Dfloat win;
  Dfloat b=4.14;
  int int_j,end_j;

  for (int j =0;j<8;j++){
    sinc_wx[j] = ZERO;
    sinc_wy[j] = ZERO;
    sinc_wz[j] = ZERO;
  }


  if (src_flag_r[is] < 0) {
    int_j = 3; end_j = 4;
  } else if (src_flag_r[is] > 0){
    int_j = 4; end_j = 3;
    }
  
  // if (freeSurf) {

  //   //b = 6.31;
  //   for (int j = -3;j<=nshift[is];j++){
  //   arg = pi * (j - 0.5);
  //   win = bessi0( b * sqrt(1 - pow((j-0.5) / 4.0,2.0))) / bessi0(b);
  //   sinc_wx[3 + j] = win * sin(arg) / (GDomain->Dx() * arg);
  //   sinc_wy[3 + j] = win * sin(arg) / (GDomain->Dy() * arg);
  //   sinc_wz[3 + j] = win * sin(arg) / (GDomain->Dz() * arg);
  // }

  //   int np = 3 + nshift[is];
  //   for (int j =nshift[is];j<=4;j++){
  //   arg = pi * (j - 0.5);
  //   win = bessi0( b * sqrt(1 - pow((j-0.5) / 4.0,2.0))) / bessi0(b);
  //   sinc_wx[np - j + nshift[is]] -= win * sin(arg) / (GDomain->Dx() * arg);
  //   sinc_wy[np - j + nshift[is]] -= win * sin(arg) / (GDomain->Dy() * arg);
  //   sinc_wz[np - j + nshift[is]] -= win * sin(arg) / (GDomain->Dz() * arg);
      
  //   }

  // } else {
    
  for (int j = -int_j;j<=end_j;j++){
    idx[int_j + j] = j;
    arg = pi * (j + 0.5 * src_flag_r[is]);
    win = bessi0( b * sqrt(1 - pow((j + 0.5 * src_flag_r[is]) / 4.0,2.0))) / bessi0(b);
    sinc_wx[int_j + j] = win * sin(arg) / (GDomain->Dx() * arg);
    sinc_wy[int_j + j] = win * sin(arg) / (GDomain->Dy() * arg);
    sinc_wz[int_j + j] = win * sin(arg) / (GDomain->Dz() * arg);
  }

  // }
  
}



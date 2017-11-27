#include "definitions.hpp"
#include "geometry3D.hpp"
#include "model.hpp"
#include "show.hpp"
#include "sdm.hpp"

int main () {

  VecF Ilim,Flim;
  VecI Nelem;
  geometry3D *Gdomain;  // Global Domain
  MODEL *model;
  Show show;
  

  Dfloat *SubMu,*SubRho,*SubLamb;

  Ilim = {0.0,0.0,0.0};
  Flim = {3100.0,3100.0,3100.0};

  Nelem = {31,31,31};
  
  Gdomain = new geometry3D(Ilim,Flim,Nelem);

  std::cout<<"DX "<<Gdomain->Dx()<<std::endl;
  std::cout<<"DY "<<Gdomain->Dy()<<std::endl;
  std::cout<<"DZ "<<Gdomain->Dz()<<std::endl;

  
  std::cout<<"L_NodX "<<Gdomain->L_NodeX()<<std::endl;
  std::cout<<"HALO_NodY "<<Gdomain->HALO_NodeY()<<std::endl;
  std::cout<<"HALO_NodZ "<<Gdomain->HALO_NodeZ()<<std::endl;
  

  VecI GNod = {Gdomain->L_NodeX(),Gdomain->L_NodeY(),Gdomain->L_NodeZ()};
  VecI SubN = {1,1,1};
  VecI SubNodes = {Gdomain->HALO_NodeX() / SubN.x, \
		   Gdomain->HALO_NodeY() / SubN.y, \
		   Gdomain->HALO_NodeZ() / SubN.z};

  std::cout<<"Nodes_Subdomain X "<<SubNodes.x<<std::endl;
  std::cout<<"Nodes_Subdomain Y "<<SubNodes.y<<std::endl;
  std::cout<<"Nodes_Subdomain Z "<<SubNodes.z<<std::endl;

  char FileVP[] = "src/example/VP.bin";
  char FileVS[] = "src/example/VS.bin";
  char FileR[] = "src/example/RHO.bin";

  model = new MODEL(FileVP,FileVS,FileR, \
		    GNod,SubNodes);


  std::cout<<"Nodes model Local "<<model->HALO_NodeX()<<std::endl;
  SubMu = new Dfloat[model->HALO_Node()];
  SubRho = new Dfloat[model->HALO_Node()];
  SubLamb = new Dfloat[model->HALO_Node()];

VecI subi = {0,0,0};
  model->SubModel(subi,SubRho,SubMu,SubLamb);

  //show.print(SubLamb,model->HALO_NodeX(),model->HALO_NodeY(),model->HALO_NodeZ());
  
  VecF GI = {Gdomain->CoorX(0),Gdomain->CoorY(0),Gdomain->CoorZ(0)};
  VecF GF = {Gdomain->CoorX(Gdomain->HALO_NodeX()-1),\
    Gdomain->CoorY(Gdomain->HALO_NodeY()-1),\
    Gdomain->CoorZ(Gdomain->HALO_NodeZ()-1)};

  Dfloat dt = 0.001;
  Dfloat f0 = 4.0;
  Dfloat t = 1.0;
  int  nt = (int) t / dt;

  SDM sdm(GI,GF,GNod,GI,GF,SubNodes,f0,dt,subi,SubN);   

  sdm.ModelRead(SubRho,"RHO");
  sdm.ModelRead(SubMu,"MU");
  sdm.ModelRead(SubLamb,"LAMB");
  int a=sdm.CFL();
  sdm.InitVar(ZERO);

  VecF post = {3100.0/2.0,3100.0/2.0,3100.0/2.0};
  VecI pos_src = Gdomain->FindNode(post);

  printf("%d\t%d\t%d\n",pos_src.x,pos_src.y,pos_src.z);

  VecI loc = sdm.SFindNode(pos_src);

  printf("%d\t%d\t%d\n",loc.x,loc.y,loc.z);

getchar();
  int kk = 0;
  for (int k = 0; k<nt; ++k){


    printf("Time step : %d of %d \n",k,nt);
  
    Dfloat source = 100000000.0 * sdm.source(1,k,0.25);
    //printf("%f\n",source);
    sdm.FDSII();
    sdm.FDSXY();
    sdm.FDSXZ();
    sdm.FDSYZ();

// Source #########
    sdm.AddVal(pos_src,"SXX",source);
    //printf("%f\n",sdm.GetVal(pos_src, "SYY"));
    sdm.AddVal(pos_src,"SYY",source);
    sdm.AddVal(pos_src,"SZZ",source);


    sdm.FDVX();
    sdm.FDVY();
    sdm.FDVZ();
    
    if (kk == k){
      sdm.print("VZ",k);
      //sdm.print("VY",k);
      //sdm.print("SXX",k);
      kk += 10;
	}

  }

 
  

return 0;
}

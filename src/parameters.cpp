#include "parameters.hpp"


PAR::PAR(int argi, char** param,char *b){

  char *nameFile = FileParam(argi,param,b);

  R.open(nameFile);

  if (R.is_open()) {

    std::cout<<"\n### PARAMETERS FILE: Read ###"<<std::endl;
  } else {

    std::cout<<"PARAMETERS FILE: Read error"<<std::endl;
  }
}

 


PAR::~PAR(){
  R.close();
}


std::string PAR::ParamReturn(std::string namePar){

  std::string line,var;

  while ( std::getline (R,line) )
      {

	if (line.find(namePar) == 0){
	  //std::cout << line << "\n"
	  int pos = line.find(" ");
	  var = line.substr(pos+1);
	  R.seekg(0, std::ios::beg);
	  break;
	}
	
      }

  return var;  
  
}


// FUNCTION TO READ FILE PARAMETERS
char * PAR::FileParam(int argi, char** param,char *b){
  char *value;
  char *c;
  
  if (argi > 0) {
    for (int i = 0; i<argi; ++i){
      if (strcmp(b,param[i]) == 0){
	value = (param[i+1]);
      }
    }
  }
  return value;
}

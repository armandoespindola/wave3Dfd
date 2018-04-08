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

#include "parameters.hpp"


PAR::PAR(int argi, char** param,char *b){

  char *nameFile = FileParam(argi,param,b);

  R.open(nameFile);

  if (R.is_open()) {

    // std::cout<<"\n### PARAMETERS FILE: Read ###"<<std::endl;
  } else {

    std::cout<<"PARAMETERS FILE: Read error"<<std::endl;
  }
}

 


PAR::~PAR(){
  R.close();
}


std::string PAR::ParamReturn(std::string namePar){

  std::string line,var;

  while ( std::getline(R,line) )
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

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

#ifndef __PARAMETERS__
#define __PARAMETERS__

#include "definitions.hpp"

// CLASS PARAMETERS
// READ THE PARAMETERS ON RUNTIME

class PAR {

protected:
  std::ifstream R;
	
public:
	
	PAR(int argi, char** param,char *b);
	~PAR();
  std::string ParamReturn(std::string namePar);
  char *FileParam(int argi, char** param,char* b);

};






#endif

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

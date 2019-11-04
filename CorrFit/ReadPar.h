#ifndef _CORRFIT_READPAR_
#define _CORRFIT_READPAR_
#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include "TString.h"

typedef TString STR;

// Ecxeption values
#define RP_Exception_UnknownException 0
#define RP_Exception_NoSuchParamter   1
#define RP_Exception_NoParFile        2

struct struct_option 
{
  STR keyword;
  STR value;
};

typedef struct struct_option option;
typedef std::vector<option> VOPT;

class ReadPar 
{
 private:
  char *fname;
  VOPT options;
  
 public:
  ReadPar(); // Default constructor
  ReadPar(const char *aFName);
  
  int readFile(const char *aFName) throw(int); 
  int printOptions();
  STR getPar(const char *name) throw(STR);
  int parseOptString(const char *aString) throw(int);
  
};

#endif

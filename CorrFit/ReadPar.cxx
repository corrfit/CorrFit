#include "CFGlobal.h"
#include "ReadPar.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iosfwd>

ReadPar::ReadPar()
{
  fname = 0;
}

ReadPar::ReadPar(const char *aFName)
{
  fname = strdup(aFName);
  readFile(aFName);
}

int ReadPar::readFile(const char *aFName) throw (int)
{
  option read_opt;
  STR buf_str;
  char buff[200];
  char dummy[200];

  std::ifstream       infile(aFName);
  std::istringstream  *instream;
  
  if (!infile.is_open())
    throw RP_Exception_NoParFile;

  instream = new std::istringstream(buff);

  while (!infile.eof())
    {
      infile.getline(buff, 200);
      instream->clear();
      instream->str(buff);

      memset(dummy,0,200);
      *instream >> dummy;
      
      PRINT_DEBUG_3("Read " << dummy);;
      read_opt.keyword = dummy;
      memset(dummy,0,200);
      *instream >> dummy;

      PRINT_DEBUG_3("Read " << dummy);
      if (strstr(dummy,"="))
	{
	  dummy[0]='\0';
	  
	  memset(dummy,0,200);
	  *instream >> dummy;
	  PRINT_DEBUG_3("Read " << dummy);
	  
	  read_opt.value = dummy;
	  options.push_back(read_opt);
	}
      
      //      delete instream;
    }
  infile.close();

  return 0;
}

int ReadPar::parseOptString(const char *aString) throw(int)
{
  option read_opt;
  STR buf_str;
  char buff[200];
  char dummy[200];
  char *cpos;
  char *cendpos;
  VOPT::iterator c;
  int iterc;

  char *instr = strdup(aString);
  cpos = instr;
  cendpos = strstr(cpos, "=");

  std::istringstream  *instream;
  instream = new std::istringstream(buff);
  
  while ((cpos) && (cendpos) && (cpos < instr+strlen(instr)))
    {
      strncpy(buff, cpos, cendpos - cpos);
      buff[cendpos  - cpos] = '\0';
      instream->clear();
      instream->str(buff);

      memset(dummy,0,200);
      *instream >> dummy;
      
      PRINT_DEBUG_3("ParString Read " << dummy);;
      read_opt.keyword = dummy;

      cpos = cendpos;
      cendpos = strstr(cpos, ",");
      if (cendpos) {
	strncpy(buff, cpos+1, cendpos - cpos - 1);
	buff[cendpos - cpos - 1] = '\0';
	cpos = cendpos+1;
	cendpos = strstr(cpos,"=");
      }
      else {
	strcpy(buff, cpos+1);
      }
      
      instream->clear();
      instream->str(buff);

      memset(dummy,0,200);
      *instream >> dummy;

      PRINT_DEBUG_3("ParString Read " << dummy);
	  
      read_opt.value = dummy;
      for (c=options.begin(), iterc=0; c != options.end(); c++, iterc++)
	if (c->keyword == read_opt.keyword)
	  {
	    PRINT_DEBUG_2("Substituting value " << read_opt.value << " for keyword " << c->keyword);
	    //	    *c->value = read_opt.value;
	    options[iterc] = read_opt;
	    break;
	  }
      if (c == options.end()) 
	options.push_back(read_opt);
    }
  

  return 0;
}


int ReadPar::printOptions()
{
  VOPT::iterator c;

  for (c=options.begin(); c != options.end(); c++)
    PRINT_DEBUG_3("Keyword: " << c->keyword << " Value: " << c->value);

  return 0;
}

STR ReadPar::getPar(const char *name) throw(STR)
{
  VOPT::iterator c;
  STR pname(name);

  for (c=options.begin(); c != options.end(); c++)
    if (c->keyword == pname)
      {
	PRINT_DEBUG_2("Returning value " << c->value << " for keyword " << c->keyword);
	return c->value;
      }
  throw *(new STR(name));
  //  throw RP_Exception_NoSuchParamter;
  return TString("");
}



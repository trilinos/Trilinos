/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <cstring>

#include <fei_sstream.hpp>

#include <test_utils/FEData.hpp>

#undef fei_file
#define fei_file "FEData.cpp"

#include <fei_ErrMacros.hpp>

int FEData::parameters(int numParams, char** params)
{
  const char* param = snl_fei::getParamValue("debugOutput",
						    numParams,params);
  if (param != NULL){
    setDebugLog(1, param);
  }

  dbgOut() << "parameters" << FEI_ENDL
	   << "   numParams: " << numParams << FEI_ENDL;
  for(int i=0; i<numParams; i++) {
    dbgOut() << "      param "<<i<<": '" << params[i] << "'" << FEI_ENDL;
  }

  return(0);
}

int FEData::setDebugLog(int debugOutputLevel, const char* path)
{
  delete [] dbgPath_;
  dbgPath_ = NULL;

  if (dbgFileOpened_ == true) return(0);

  if (path != NULL) {
    dbgPath_ = new char[strlen(path)+1];
    std::strcpy(dbgPath_, path);
  }
  else {
    dbgPath_ = new char[2];
    std::strcpy(dbgPath_, ".");
  }

  debugOutputLevel_ = debugOutputLevel;

  if (debugOutputLevel_ <= 0) {
    dbgOStreamPtr_ = NULL;
  }
  else {
    if (dbgOStreamPtr_ != NULL) delete dbgOStreamPtr_;
    dbgOStreamPtr_ = NULL;

    FEI_OSTRINGSTREAM fname;
    fname << dbgPath_<<"/FEData."<<numProcs_<<"."<<localProc_;
    dbgFStreamPtr_ = new FEI_OFSTREAM(fname.str().c_str());
    dbgFileOpened_ = true;
    dbgOStreamPtr_ = dbgFStreamPtr_;
  }

  return(0);
}


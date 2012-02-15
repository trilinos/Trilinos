/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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


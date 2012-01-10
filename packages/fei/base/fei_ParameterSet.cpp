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


#include <fei_macros.hpp>
#include <fei_ParameterSet.hpp>

fei::ParameterSet::ParameterSet()
  : params_(NULL)
{
  params_ = new std::vector<const Param*>;
}

fei::ParameterSet::~ParameterSet()
{
  const_iterator iter = begin(), iter_end = end();
  for(; iter != iter_end; ++iter) {
    delete &(*iter);
  }

  delete params_; params_ = NULL;
}

int fei::ParameterSet::getIntParamValue(const char* name,
					int& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::INT) return(-1);

  paramValue = param->getIntValue();
  return(0);
}

int fei::ParameterSet::getDoubleParamValue(const char* name,
					   double& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() == fei::Param::INT) {
    paramValue = param->getIntValue();
  }
  else if (param->getType() == fei::Param::DOUBLE) {
    paramValue = param->getDoubleValue();
  }
  else return(-1);

  return(0);
}

int fei::ParameterSet::getStringParamValue(const char* name,
					   std::string& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::STRING) return(-1);

  paramValue = param->getStringValue();
  return(0);
}

int fei::ParameterSet::getBoolParamValue(const char* name,
					 bool& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::BOOL) return(-1);

  paramValue = param->getBoolValue();
  return(0);
}

int fei::ParameterSet::getVoidParamValue(const char* name,
					 const void*& paramValue) const
{
  const fei::Param* param = get(name);
  if (param == NULL) return(-1);

  if (param->getType() != fei::Param::VOID) return(-1);

  paramValue = param->getVoidValue();
  return(0);
}

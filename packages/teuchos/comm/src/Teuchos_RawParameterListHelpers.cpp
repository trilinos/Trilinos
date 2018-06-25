// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "Teuchos_RawParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommHelpers.hpp"

TEUCHOSCOMM_LIB_DLL_EXPORT
void Teuchos::updateParametersAndBroadcast(
  const Ptr<ParameterList> &inParamList,                                     
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm,
  int root,
  bool overwrite)
{
  
  if (comm.getSize()==1) {
    if(overwrite) paramList->setParameters(*inParamList);
    else paramList->setParametersNotAlreadySet(*inParamList);
  }
  else {
    if (comm.getRank()==root) {
      XMLParameterListWriter w;
      std::string xmlString = toString(w.toXML(*inParamList));
      int strsize = static_cast<int>(xmlString.size());
      broadcast<int, int>(comm, root, &strsize);
      broadcast<int, char>(comm, root, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
    else {
      int strsize;
      broadcast<int, int>(comm, root, &strsize);
      std::string xmlString;
      xmlString.resize(strsize);
      broadcast<int, char>(comm, root, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
  }
}

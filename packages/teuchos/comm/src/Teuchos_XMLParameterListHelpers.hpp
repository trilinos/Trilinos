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

#ifndef TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP
#define TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP


/*! \file Teuchos_XMLParameterListHelpers.hpp \brief Additional ParameterList
  XML helper functions including parallel support.
*/


#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_Comm.hpp"


namespace Teuchos {


/** \brief On processor rank = 0, reads XML parameters from a file 
 * and broadcasts them to all other processors. Then updates the 
 * given parameter list with these values.
 *
 * \param xmlFileName [in] The file name containing XML parameter list
 * specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>xmlFileName</tt> will be set or overide those in
 * <tt>*paramList</tt>.
 *
 * \param comm [in] A Comm object used to broadcast the xml.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT
void updateParametersFromXmlFileAndBroadcast(
  const std::string &xmlFileName,
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm
  );


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP

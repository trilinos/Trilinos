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

#ifndef TEUCHOS_YAML_PARAMETER_LIST_CORE_HELPERS_HPP
#define TEUCHOS_YAML_PARAMETER_LIST_CORE_HELPERS_HPP


/*! \file Teuchos_YamlParameterListCoreHelpers.hpp \brief Simple helper functions
     that make it easy to read and write Yaml to and from a parameterlist.
*/


#include "Teuchos_ParameterList.hpp"


namespace Teuchos {


/** \brief Reads Yaml parameters from a file and updates those already in the
 * given parameter list.
 *
 * \param yamlFileName [in] The file name containing Yaml parameter list
 * specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>yamlFileName</tt> will be set or overide those in
 * <tt>*paramList</tt>.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void updateParametersFromYamlFile(
  const std::string &yamlFileName,
  const Ptr<ParameterList> &paramList
  );

/** \brief Reads Yaml parameters from a file and return them in a new parameter
 * list.
 *
 * \param yamlFileName [in] The file name containing Yaml parameter list
 * specification.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromYamlFile(const std::string &yamlFileName);

/** \brief Reads Yaml parameters from a std::string and updates those already in the
 * given parameter list.
 *
 * \param yamlStr [in] String containing Yaml parameter list specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>yamlStr</tt> will be set or override (or not) those in
 * <tt>*paramList</tt> depending on the <tt>overwrite</tt> parameter.
 *
 * \param overwrite [in] If true, parameters and sublists in the <tt>yamlStr</tt> 
 * will override those in <tt>paramList</tt>.  If false, any value set in 
 * <tt>paramList</tt> will be kept, only values not set will be updated.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void updateParametersFromYamlString(
  const std::string &yamlStr,
  const Ptr<ParameterList> &paramList, 
  bool overwrite,
  const std::string& name = ""
  );

TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT 
void updateParametersFromYamlCString(
  const char* const data,
  const Teuchos::Ptr<Teuchos::ParameterList>& paramList,
  bool overwrite
  );



/** \brief Reads Yaml parameters from a std::string and return them in a new
 * parameter list.
 *
 * \param yamlStr [in] String containing Yaml parameter list specification.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromYamlString(const std::string &yamlStr);


/** \brief Write parameters and sublists in Yaml format to an std::ostream.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written to file.
 *
 * \param yamlOut [in] The stream that will get the Yaml output.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void writeParameterListToYamlOStream(
  const ParameterList &paramList,
  std::ostream &yamlOut
  );


/** \brief Write parameters and sublist to an Yaml file.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written to file.
 *
 * \param yamlFileName [in] The file name that will be create to contain the
 * Yaml version of the parameter list specification.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void writeParameterListToYamlFile(
  const ParameterList &paramList,
  const std::string &yamlFileName
  );

} // namespace Teuchos


#endif // TEUCHOS_Yaml_PARAMETER_LIST_CORE_HELPERS_HPP

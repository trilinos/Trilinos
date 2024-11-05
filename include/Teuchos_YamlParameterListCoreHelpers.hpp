// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_XML_PARAMETER_LIST_CORE_HELPERS_HPP
#define TEUCHOS_XML_PARAMETER_LIST_CORE_HELPERS_HPP


/*! \file Teuchos_XMLParameterListCoreHelpers.hpp \brief Simple helper functions
     that make it easy to read and write XML to and from a parameterlist.
*/


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DependencySheet.hpp"


namespace Teuchos {


/** \brief Reads XML parameters from a file and updates those already in the
 * given parameter list.
 *
 * \param xmlFileName [in] The file name containing XML parameter list
 * specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>xmlFileName</tt> will be set or overide those in
 * <tt>*paramList</tt>.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void updateParametersFromXmlFile(
  const std::string &xmlFileName,
  const Ptr<ParameterList> &paramList
  );


/** \brief Reads XML parameters from a file and return them in a new parameter
 * list.
 *
 * \param xmlFileName [in] The file name containing XML parameter list
 * specification.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlFile(const std::string &xmlFileName);


/** \brief Reads XML parameters from a file and return them in a new parameter
 * list.
 *
 * \param xmlFileName [in] The file name containing XML parameter list
 * specification.
 *
 * \param depSheet [out] The Dependency Sheet into which Dependencies should be
 * placed.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlFile(const std::string &xmlFileName,
  RCP<DependencySheet> depSheet);


/** \brief Reads XML parameters from a std::string and updates those already in the
 * given parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>xmlStr</tt> will be set or override (or not) those in
 * <tt>*paramList</tt> depending on the <tt>overwrite</tt> parameter.
 *
 * \param overwrite [in] If true, parameters and sublists in the <tt>xmlStr</tt> 
 * will override those in <tt>paramList</tt>.  If false, any value set in 
 * <tt>paramList</tt> will be kept, only values not set will be updated.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void updateParametersFromXmlString(
  const std::string &xmlStr,
  const Ptr<ParameterList> &paramList, 
  bool overwrite = true
  );


/** \brief Reads XML parameters from a std::string and return them in a new
 * parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlString(const std::string &xmlStr);


/** \brief Reads XML parameters from a std::string and return them in a new
 * parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 * \param depSheet [in] The Dependency Sheet into which Dependencies should be
 * placed.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlString( const std::string &xmlStr,
  RCP<DependencySheet> depSheet);


/** \brief Write parameters and sublists in XML format to an std::ostream.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written to file.
 *
 * \param xmlOut [in] The stream that will get the XML output.
 *
 * \param depSheet [in] The Dependency Sheet which should be written out.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void writeParameterListToXmlOStream(
  const ParameterList &paramList,
  std::ostream &xmlOut,
  RCP<const DependencySheet> depSheet = null
  );


/** \brief Write parameters and sublist to an XML file.
 *
 * \param paramList [in] Contains the parameters and sublists that will be
 * written to file.
 *
 * \param xmlFileName [in] The file name that will be create to contain the
 * XML version of the parameter list specification.
 *
 * \param depSheet [in] The Dependency Sheet which should be written out.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT
void writeParameterListToXmlFile(
  const ParameterList &paramList,
  const std::string &xmlFileName,
  RCP<const DependencySheet> depSheet=null
  );


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_CORE_HELPERS_HPP

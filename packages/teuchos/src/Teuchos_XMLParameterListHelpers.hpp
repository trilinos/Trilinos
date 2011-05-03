// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP
#define TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP


/*! \file Teuchos_XMLParameterListHelpers.hpp \brief Simple helper functions
     that make it easy to read and write XML to and from a parameterlist.
*/


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_Comm.hpp"


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
TEUCHOS_LIB_DLL_EXPORT void updateParametersFromXmlFile(
  const std::string &xmlFileName,
  ParameterList *paramList
  );

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
  ParameterList *paramList,
  const Comm<int> &comm
  );

/** \brief Reads XML parameters from a file and return them in a new parameter 
 * list.
 *
 * \param xmlFileName [in] The file name containing XML parameter list
 * specification.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT 
RCP<ParameterList> getParametersFromXmlFile( const std::string &xmlFileName );

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
TEUCHOS_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlFile(const std::string &xmlFileName,
  RCP<DependencySheet> depSheet);


/** \brief Reads XML parameters from a std::string and updates those already in the
 * given parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from the file <tt>xmlStr</tt> will be set or overide those in
 * <tt>*paramList</tt>.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT
void updateParametersFromXmlString(
  const std::string &xmlStr,
  ParameterList *paramList
  );


/** \brief Reads XML parameters from a std::string and return them in a new
 * parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT
RCP<ParameterList> getParametersFromXmlString( const std::string &xmlStr );

/** \brief Reads XML parameters from a std::string and return them in a new
 * parameter list.
 *
 * \param xmlStr [in] String containing XML parameter list specification.
 * \param depSheet [in] The Dependency Sheet into which Dependencies should be
 * placed.
 * 
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT
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
TEUCHOS_LIB_DLL_EXPORT
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
TEUCHOS_LIB_DLL_EXPORT
void writeParameterListToXmlFile(
  const ParameterList &paramList,
  const std::string &xmlFileName,
  RCP<const DependencySheet> depSheet=null
  );


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP

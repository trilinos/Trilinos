// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
 * from the file <tt>xmlFileName</tt> will be set or override (or not) those in
 * <tt>*paramList</tt> depending on the <tt>overwrite</tt> parameter.
 *
 * \param comm [in] A Comm object used to broadcast the xml.
 *
 * \param overwrite [in] If true, parameters and sublists in the <tt>xmlStr</tt> 
 * will override those in <tt>paramList</tt>.  If false, any value set in 
 * <tt>paramList</tt> will be kept, only values not set will be updated.
 *
 * \relates ParameterList
 */
TEUCHOSCOMM_LIB_DLL_EXPORT
void updateParametersFromXmlFileAndBroadcast(
  const std::string &xmlFileName,
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm,
  bool overwrite = true
  );


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP

// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_RAW_PARAMETER_LIST_HELPERS_HPP
#define TEUCHOS_RAW_PARAMETER_LIST_HELPERS_HPP


/*! \file Teuchos_RAWParameterListHelpers.hpp \brief Additional ParameterList
  helper functions including parallel support.
*/
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Comm.hpp"


namespace Teuchos {


/** \brief On processor rank = root, broadcast the inParamList
 * to all other processors. Then update the given parameter list with 
 * these values.
 *
 * \param inParamList [in] On input, <tt>*inParamList</tt> may be empty or
 * contain some parameters and sublists. This list is only used on rank = root.
 *
 * \param paramList [in/out] On input, <tt>*paramList</tt> may be empty or
 * contain some parameters and sublists. On output, parameters and sublist
 * from rank = root will be set or override (or not) those in
 * <tt>*paramList</tt> depending on the <tt>overwrite</tt> parameter.
 *
 * \param comm [in] A Comm object used to broadcast the xml.
 *
 * \param root [in] The rank which will broadcast the list.
 *
 * \param overwrite [in] If true, parameters and sublists in the <tt>inParamList</tt> 
 * will override those in <tt>paramList</tt>.  If false, any value set in 
 * <tt>paramList</tt> will be kept, only values not set will be updated.
 *
 * \relates ParameterList
 */
TEUCHOSCOMM_LIB_DLL_EXPORT
void updateParametersAndBroadcast(
  const Ptr<ParameterList> &inParamList,                                     
  const Ptr<ParameterList> &ParamList,
  const Comm<int> &comm,
  int root,
  bool overwrite = true
  );


} // namespace Teuchos


#endif // TEUCHOS_XML_PARAMETER_LIST_HELPERS_HPP

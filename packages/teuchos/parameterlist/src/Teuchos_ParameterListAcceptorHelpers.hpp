// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETER_LIST_ACCEPTOR_HELPERS_HPP
#define TEUCHOS_PARAMETER_LIST_ACCEPTOR_HELPERS_HPP


#include "Teuchos_ParameterListAcceptor.hpp"


namespace Teuchos {


class ParameterListAcceptor;


/** \brief Pretty print the valid parameters from a ParameterListAccpetor
 * object.
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void printValidParameters( const ParameterListAcceptor &paramListAccpetor,
  std::ostream &out, const bool showDoc = true );


} // end namespace Teuchos


#endif // TEUCHOS_PARAMETER_LIST_ACCEPTOR_HELPERS_HPP

// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterListAcceptorHelpers.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_ParameterList.hpp"


void Teuchos::printValidParameters( const ParameterListAcceptor &paramListAccpetor,
  std::ostream &out, const bool showDoc )
{
  typedef ParameterList::PrintOptions PLPrintOptions;
  paramListAccpetor.getValidParameters()->print(
    out, PLPrintOptions().indent(2).showTypes(true).showDoc(showDoc)
    );
}

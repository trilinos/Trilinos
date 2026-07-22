// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterListNonAcceptor.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teuchos {


// Overridden from ParameterListAcceptor


void ParameterListNonAcceptor::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*this->getValidParameters());
  setMyParamList(paramList);
}


RCP<const ParameterList>
ParameterListNonAcceptor::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    validPL = parameterList();
  }
  return validPL;
}


} // end namespace Teuchos

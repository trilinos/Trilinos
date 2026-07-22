// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Teuchos {


// Overridden from ParameterListAcceptor


RCP<ParameterList>
ParameterListAcceptorDefaultBase::getNonconstParameterList()
{
  return paramList_;
}


RCP<ParameterList>
ParameterListAcceptorDefaultBase::unsetParameterList()
{
  RCP<ParameterList> tempParamList = paramList_;
  paramList_ = Teuchos::null;
  return tempParamList;
}


RCP<const ParameterList>
ParameterListAcceptorDefaultBase::getParameterList() const
{
  return paramList_;
}


} // namespace Teuchos

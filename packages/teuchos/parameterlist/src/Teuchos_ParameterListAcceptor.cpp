// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DependencySheet.hpp"


namespace Teuchos {


ParameterListAcceptor::~ParameterListAcceptor()
{}


Teuchos::RCP<const Teuchos::ParameterList>
ParameterListAcceptor::getParameterList() const
{
  return const_cast<ParameterListAcceptor*>(this)->getNonconstParameterList();
}


Teuchos::RCP<const Teuchos::ParameterList>
ParameterListAcceptor::getValidParameters() const
{
  return Teuchos::null;
}

RCP<const DependencySheet>
ParameterListAcceptor::getDependencies() const
{
  return null;
}


} // end namespace Teuchos

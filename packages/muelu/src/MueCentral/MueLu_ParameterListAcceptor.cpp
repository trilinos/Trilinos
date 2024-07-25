// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include "MueLu_ParameterListAcceptor.hpp"

// TODO See also: Teuchos::ParameterListAcceptor, Teko::Clonable

namespace MueLu {

void printParameterListOptions(std::ostream& os, const Teuchos::ParameterList& p) {
  p.print(os, Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));
  os << std::endl;
}

ParameterListAcceptor::ParameterListAcceptor() {}

ParameterListAcceptorImpl::ParameterListAcceptorImpl() {}

void ParameterListAcceptorImpl::SetParameterList(const Teuchos::ParameterList& paramList) {
  // This call is only for cosmetic reasons.
  // If one calls SetParameterList before GetParameterList, that would mean
  // that paramList_ has not been initialized yet. Therefore, the parameter
  // would be put in it in the order user provided, and not in the order of
  // valid parameter list. We'd like to have consistency in the order, so
  // we do this extra call, which is no-op if paramList_ has already been
  // initialized.
  paramList_ = GetParameterList();

  paramList_.setParameters(paramList);

  // Validate and add defaults parameters.
  Teuchos::RCP<const Teuchos::ParameterList> validParamList = GetValidParameterList();
  if (validParamList != Teuchos::null) {
    paramList_.validateParametersAndSetDefaults(*validParamList);
  } else {
    // Teuchos::null means that GetValidParameterList() not implemented.
    // As such, we skip validation and have not way to set default values,
    // which is potentially dangerous
  }
}

const Teuchos::ParameterList& ParameterListAcceptorImpl::GetParameterList() const {
  // The returned list always has an entry for each valid parameter.
  // Therefore, there is not need to test if a parameter is present before getting it.
  if (paramList_.numParams() == 0) {
    // Set paramList_ to the default list
    Teuchos::RCP<const Teuchos::ParameterList> validParamList = GetValidParameterList();
    if (validParamList != Teuchos::null) {
      // Instead of simply doing
      //   paramList_ = *validParamList;
      // we use more complicated Teuchos calls, because we would like to
      // have [default] values in the beginning
      paramList_.validateParametersAndSetDefaults(*validParamList);
    }

  } else {
    // We are sure that the list has all the valid parameters defined
    // because the parameter list validation process adds the default
    // values to the user list
  }

  return paramList_;
}

void ParameterListAcceptorImpl::SetParameter(const std::string& name, const ParameterEntry& entry) {
  Teuchos::ParameterList paramList;
  paramList.setEntry(name, entry);
  SetParameterList(paramList);  // This forces revalidation of the list
}

const ParameterEntry& ParameterListAcceptorImpl::GetParameter(const std::string& name) const {
  return GetParameterList().getEntry(name);
}

void ParameterListAcceptorImpl::GetDocumentation(std::ostream& os) const {
  // default implementation

  Teuchos::RCP<const Teuchos::ParameterList> validParamList = GetValidParameterList();
  if (validParamList == Teuchos::null) {
    os << "## Documentation not available:" << std::endl;
    return;
  }

  os << "## Parameters:" << std::endl;
  printParameterListOptions(os, *validParamList);

  os << "## Fully described default method:" << std::endl;
  validParamList->print(os, 2, true, false);
  os << std::endl;
}

}  // namespace MueLu

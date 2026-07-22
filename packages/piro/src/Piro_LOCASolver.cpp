// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_LOCASolver_Def.hpp"

namespace Piro {

// Explicit template instantiations
// LOCA currently only supports Scalar = double

template class LOCASolver<double>;

template Teuchos::RCP<LOCASolver<double> > observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> > &adjointModel,
    const Teuchos::RCP<Piro::ObserverBase<double> > &observer);

} // namespace Piro


#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_toString.hpp"

Piro::Detail::ModelEvaluatorParamName::ModelEvaluatorParamName(
    const Teuchos::RCP<const Teuchos::Array<std::string> > &p_names) :
  p_names_(p_names)
{
  if (Teuchos::is_null(p_names_)) {
    type_ = Default;
  } else if (p_names_->size() == Teuchos::OrdinalTraits<Teuchos_Ordinal>::one()) {
    type_ = OneShared;
  } else {
    type_ = FullList;
  }
}

std::string
Piro::Detail::ModelEvaluatorParamName::operator()(Teuchos_Ordinal k) const
{
  switch (type_) {
    case Default:
      return "Parameter " + Teuchos::toString(k);
    case OneShared:
      return p_names_->front();
    case FullList:
      return (*p_names_)[k];
  }

  TEUCHOS_TEST_FOR_EXCEPT(true);
  return std::string();
}

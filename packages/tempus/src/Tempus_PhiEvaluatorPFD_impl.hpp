//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorPFD_impl_hpp
#define Tempus_PhiEvaluatorPFD_impl_hpp

#include "Tempus_PhiEvaluatorPFD.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluatorPFD_decl.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorPFD<Scalar>::getValidParameters() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  //pl->set(
  //    "var", default,
  //    "'var' sets the var.  "
  //    "'opt1' - will do this.  "
  //    "'opt2' - will do that!");

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
void PhiEvaluatorPFD<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>>,
                                         int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> rhs_b)
{
  //TODO
}
  
// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorPFD<Scalar> > createPhiEvaluatorPFD(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto phi = rcp(new PhiEvaluatorPFD<Scalar>());
  phi->setName("From createPhiEvaluatorPFD");

  if (pl == Teuchos::null || pl->numParams() == 0) return phi;

  pl->validateParametersAndSetDefaults(*phi->getValidParameters());

  phi->setName(pl->name());
  //phi->setThing(pl->get("Thing", "default"));

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorPFD_impl_hpp

//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_impl_hpp
#define Tempus_PhiEvaluatorLeja_impl_hpp

#include "Tempus_PhiEvaluatorLeja.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluatorLeja_decl.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorLeja<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "Leja",
      "Method to approximate the phi-function evaluation.");

  pl->set<int>(
      "Max Leja Expansion Order", 500,
      "Maximal order of the Leja polynomial used.\n"
      "\n"
      "The default is 500.");

  return pl;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorLeja<Scalar>::computeLinOpPhi(const int phi_order,
								     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
								     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order != 0,
      std::invalid_argument,
      "LinOpPhi: phi_order must be zero.");

  // TODO: Implement
  return Thyra::SolveStatus<Scalar>();
}
  
template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  maxLejaOrder_ = pl->get<int>("Max Leja Order", 500);

  // TODO: has to be set to true, only matrix exponential is implemented
  this->useAtildeForSingleRHS_ = true;

  std::cout << "\nuseAtildeForSingleRHS_: " << this->useAtildeForSingleRHS_ << std::endl;
  std::cout << "Parameter List: " << *pl << std::endl;
  std::cout << "Leja Order is " << maxLejaOrder_ << std::endl;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<PhiEvaluatorLeja<Scalar>> phi = Teuchos::rcp(new PhiEvaluatorLeja<Scalar>());
  phi->setName("From createPhiEvaluatorLeja");

  if (pl != Teuchos::null)
    phi->setPhiEvaluatorValues(pl);

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorLeja_impl_hpp

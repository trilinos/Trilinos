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

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorPFD<Scalar>::getValidParameters() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  pl->set(
      "PFD method", "CN",
      "'PDF method' determines the partial fraction decomposition used to approximate the exponential.  "
      "'IE' - uses an implicit Euler approximation (order 1).  "
      "'CN' - uses a Crank-Nicolson approximation (order 2).");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
void PhiEvaluatorPFD<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs)
{
  inArgs_lin_ = Teuchos::rcpFromRef(inArgs);
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorPFD<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
							       int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> rhs_b)
{
  // TODO: right now, hard-codes 'CN' method and k == 1. Generalize.
  
  const Scalar alpha = Scalar(1.0);
  const Scalar beta  = Scalar(0.5) * cdt;

  std::cout << "computing Phi!" << std::endl;
  
  Thyra::SolveStatus<Scalar> sStatus = this->phiLinSolv_->solveMpJ(*inArgs_lin_, phiv, rhs_b, alpha, beta);

  return sStatus;
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

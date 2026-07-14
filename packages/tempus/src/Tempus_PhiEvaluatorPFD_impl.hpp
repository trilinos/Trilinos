//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
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
      "PFD Method", "CN",
      "'PDF Method' determines the partial fraction decomposition used to approximate the exponential.  "
      "'IE' - uses an implicit Euler approximation (order 1).  "
      "'CN' - uses a Crank-Nicolson approximation (order 2).");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorPFD<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
							       const int k, const Scalar cdt,
							       const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b)
{
  // TODO: right now, hard-codes 'CN' method and k == 1. Generalize.

  TEUCHOS_TEST_FOR_EXCEPTION(
      k != 1,
      std::invalid_argument,
      "PhiEvaluatorPFD<Scalar>::computePhi is only implemented for k=1");

  // check that model and PhiLinearSolver have been initialized properly
  this->checkInitialized();

  const Scalar alpha = Scalar(1.0);
  const Scalar beta  = Scalar(0.5) * cdt;

  // TODO: we call solveMpJ without inArgs, it will use the existing Mass and Jacobian.
  //       check if there is aperformance penalty
  //Thyra::SolveStatus<Scalar> sStatus = this->phiLinSolv_->assembleAndsolveMpJ(this->inArgs_lin_, phiv, Mrhs_b, alpha, beta);
  Thyra::SolveStatus<Scalar> sStatus = this->phiLinSolv_->solveMpJ(phiv, Mrhs_b, alpha, beta);

  return sStatus;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorPFD<Scalar>::computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                                const Scalar cdt,
                                                                const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B)
{
  bool not_phi1 = (Mrhs_B.size() != 2) || (Mrhs_B[0] != Teuchos::null);

  TEUCHOS_TEST_FOR_EXCEPTION(
      not_phi1,
      std::invalid_argument,
      "PhiEvaluatorPFD<Scalar>::computePhis is only implemented for k=1");

  return this->computePhi(x, 1, cdt, Mrhs_B[1]);
}

template <class Scalar>
void PhiEvaluatorPFD<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  std::string pfdMethod = pl->get<std::string>("PFD Method", "CN");
  if (pfdMethod != "CN")
  {
    Teuchos::RCP<Teuchos::FancyOStream> l_out = this->getOStream();
    l_out->setOutputToRootOnly(0);
    *l_out << "PFD Method '" << pfdMethod << "'\n"
           << "is not yet implemented, continuing with CN!\n";
  }

  if (this->lumpMassMatrix_ == true)
  {
    Teuchos::RCP<Teuchos::FancyOStream> l_out = this->getOStream();
    l_out->setOutputToRootOnly(0);
    *l_out << "Option: 'Lump Mass Matrix' is not supported for PFD Solvers, continuing with full matrix.\n";
    this->setLumpMassMatrix(false);
  }

  std::cout << "\nParameter List: " << *pl << std::endl;
}


// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorPFD<Scalar>> createPhiEvaluatorPFD(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto phi = rcp(new PhiEvaluatorPFD<Scalar>());
  phi->setName("From createPhiEvaluatorPFD");

  if (pl != Teuchos::null)
    phi->setPhiEvaluatorValues(pl);

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorPFD_impl_hpp

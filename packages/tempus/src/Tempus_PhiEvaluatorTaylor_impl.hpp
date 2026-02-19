//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_impl_hpp
#define Tempus_PhiEvaluatorTaylor_impl_hpp

#include "Tempus_PhiEvaluatorTaylor.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorTaylor<Scalar>::getValidParameters() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set<int>(
      "Taylor Expansion Order", 4,
      "The order of the Taylor expansion used in the EPI stepper.\n"
      "\n"
      "The default is 2.");

  pl->set<bool>(
  "Zero Initial Guess", false,
  "Indicates whether to use a zero initial guess for the nonlinear\n"
  "solver when computing phi-function evaluations.");

  pl->set<std::string>("Solver Name", "Demo Solver", "Solver name for PhiEvaluator.");
  pl->set<std::string>("Predictor Stepper Type", "None", "Solver name for PhiEvaluator.");

  pl->set<bool>("Lump Mass Matrix", true, "Whether to lump the mass matrix in PhiEvaluator.");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs)
{
  inArgs_lin_ = Teuchos::rcpFromRef(inArgs);
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorTaylor<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
							       int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mrhs_b)
{
  // phi->setLumpMassMatrix(useLumpedMass_);
  this->phiLinSolv_->setLumpMassMatrix(useLumpedMass_);
  this->phiLinSolv_->computeMassMatrix(*inArgs_lin_);
  this->phiLinSolv_->computeJacobian(*inArgs_lin_);
  this->phiLinSolv_->buildK(k);
  Teuchos::RCP<Thyra::VectorBase<Scalar>> rhs_b = Mrhs_b->clone_v();
  this->phiLinSolv_->solveMass(rhs_b.ptr(), Mrhs_b);

  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // rhs_b->describe(*out, Teuchos::VERB_EXTREME);
  // Mrhs_b->describe(*out, Teuchos::VERB_EXTREME);
  //   auto vec = rhs_b;
  // auto space = vec->space();
  // int n = space->dim();

  // for (int i = 0; i < n; ++i) {
  //    std::cout << "rhs[" << i << "] = "
  //              << Thyra::get_ele(*rhs_b, i) << std::endl;
  //              std::cout << "Mrhs[" << i << "] = "
  //              << Thyra::get_ele(*Mrhs_b, i) << std::endl;
  // }

  this->phiLinSolv_->buildb(k, rhs_b);
  Atilde_ = this->phiLinSolv_->buildATilde(cdt);
  v_ = this->phiLinSolv_->buildv(Atilde_->domain());
  auto vec = this->matrixExponential(taylorExpOrder_);

  // Get the first block of the multi-vector calculated from 2x2 multi-matrix
  auto pv = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar>>(vec, true);
  auto v0 = pv->getVectorBlock(0);  // V block
  Thyra::copy(*v0, phiv.ptr());

  Thyra::SolveStatus<Scalar> sStatus;
  sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
  return sStatus;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>> PhiEvaluatorTaylor<Scalar>::matrixExponential(const int expansionOrder)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      expansionOrder < 0,
      std::invalid_argument,
      "matrixExponential: expansionOrder must be nonnegative");

  // exp(A) * v is in range(A)
  const auto rangeSpace = Atilde_->range();

  // Create tmp vector to hold result
  auto matExpTemp = Thyra::createMember(rangeSpace);

  Thyra::assign(matExpTemp.ptr(), Scalar(0));

  // Iteration vector d_0 = v
  Teuchos::RCP<Thyra::VectorBase<Scalar>> d_k = Thyra::createMember(rangeSpace);
  Thyra::assign(d_k.ptr(), *v_);

  // matExpTemp += d_0 / 0!
  Thyra::Vp_V(matExpTemp.ptr(), *d_k);

  Teuchos::RCP<Thyra::VectorBase<Scalar>> next = Thyra::createMember(rangeSpace);
  Scalar err_est;
  Scalar overflow = 0.;
  int k;
  // Iteratively compute d_k = (A^k v) / k! and add to result
  for (k = 1; k <= expansionOrder; ++k)
  {
    // next <- A * d_k
    // TODO: do we need the temp vector?
    Thyra::apply(*Atilde_, Thyra::NOTRANS, *d_k, next.ptr());

    // multiply the update by 1/k and store in d_k
    Thyra::V_StV(d_k.ptr(), Scalar(1.) / Scalar(k), *next);

    // add d_k to the final result
    Thyra::Vp_V(matExpTemp.ptr(), *d_k);

    err_est = Thyra::norm_inf(*d_k);

    if (err_est < 1e-20)
      break;
  }
  std::cout << "Norm of final update in iteration " << k << " is " << err_est << std::endl;

  matExp_v_ = matExpTemp; // This is required to wrap multivector as linearop

  return matExp_v_;
}

template<class Scalar>
void PhiEvaluatorTaylor<Scalar>::setLumpMassMatrix(bool lump)
{
  if (this->phiLinSolv_ != Teuchos::null)
  {
    std::cout << "Setting lump mass matrix to " << lump << std::endl;
    this->phiLinSolv_->setLumpMassMatrix(lump);
  }
    useLumpedMass_ = lump;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar> > createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto phi = Teuchos::rcp(new PhiEvaluatorTaylor<Scalar>());
  phi->setName("From createPhiEvaluatorTaylor");
  if (pl == Teuchos::null || pl->numParams() == 0) return phi;

  pl->validateParametersAndSetDefaults(*phi->getValidParameters());

  phi->setTaylorExpansionOrder(pl->get<int>("Taylor Expansion Order", 10));
  std::cout << "Parameter List: " << *pl << std::endl;

  // phi->setLumpMassMatrix(false);
  phi->setLumpMassMatrix(pl->get<bool>("Lump Mass Matrix", true));
  phi->setName(pl->name());
  //phi->setThing(pl->get("Thing", "default"));

  std::cout << "Taylor Expansion Order is " << phi->getTaylorExpansionOrder() << std::endl;

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorTaylor_impl_hpp

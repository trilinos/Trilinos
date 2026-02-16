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
  // return vecV;
  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // vec->describe(*out, Teuchos::VERB_EXTREME);
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

  // Identity * v = v
  Teuchos::RCP<Thyra::VectorBase<Scalar>> term = Thyra::createMember(rangeSpace);
  Thyra::assign(term.ptr(), *v_);

  // matExpTemp += term / 0!
  Thyra::Vp_V(matExpTemp.ptr(), *term);

  // Iteratively compute term = A * term (A^k v) and accumulate term/k!
  Scalar invFact = Scalar(1); // 1/k! updated each step
  for (int k = 1; k <= expansionOrder; ++k)
  {
    // term <- A * term
    Teuchos::RCP<Thyra::VectorBase<Scalar>> next = Thyra::createMember(rangeSpace);
    Thyra::apply(*Atilde_, Thyra::NOTRANS, *term, next.ptr());
    term = next;

    invFact /= Scalar(k);

    // multiply with inverse factorial
    Thyra::Vp_StV(matExpTemp.ptr(), invFact, *term);
  }

  matExp_v_ = matExpTemp; // This is required to wrap multivector as linearop

  return matExp_v_;
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

  auto test = pl->get<int>("Taylor Expansion Order", 4);
  phi->setTaylorExpansionOrder(test);
  phi->setName(pl->name());
  //phi->setThing(pl->get("Thing", "default"));

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorTaylor_impl_hpp

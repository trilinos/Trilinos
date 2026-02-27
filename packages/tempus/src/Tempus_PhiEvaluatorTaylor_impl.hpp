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
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "Taylor",
      "Method to approximate the phi-function evaluation.");

  pl->set<int>(
      "Taylor Expansion Order", 10,
      "The order of the Taylor expansion used in the EPI stepper.\n"
      "\n"
      "The default is 10.");

  return pl;
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs)
{
  inArgs_lin_ = Teuchos::rcpFromRef(inArgs);
}

template <class Scalar>
Thyra::SolveStatus<Scalar>
PhiEvaluatorTaylor<Scalar>::computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
					Scalar cdt,
					const std::vector<Teuchos::RCP<const Thyra::VectorBase<Scalar>>> Mrhs_B)
{
  int p = Mrhs_B.size() - 1;

  TEUCHOS_TEST_FOR_EXCEPTION(
      p < 1,
      std::invalid_argument,
      "computePhis: list of rhs must have at least two entries.");
  // TODO: Support p = 0? It can be done by calling computePhi with k=0, but that requires a dedicated impl.

  this->phiLinSolv_->setLumpMassMatrix(this->lumpMassMatrix_);
  this->phiLinSolv_->computeMassMatrix(*inArgs_lin_);
  this->phiLinSolv_->computeJacobian(*inArgs_lin_);

  std::vector<Teuchos::RCP<const Thyra::VectorBase<Scalar>>> rhs_B(p+1);

  // Invert the mass matrix out of the right hand sides
  // TODO: This might be more efficient to do on the combined MultiVector that will be assembled in buildb
  //       However, if Mrhs_B is sparse, it may not.
  for (int ii = 0; ii < p+1; ii++)
  {
    if (Mrhs_B[ii] != Teuchos::null)
    {
      auto Mrhs_b = Mrhs_B[ii];
      Teuchos::RCP<Thyra::VectorBase<Scalar>> rhs_b = Mrhs_b->clone_v();
      Thyra::assign(rhs_b.ptr(), 0.0);
      this->phiLinSolv_->solveMass(rhs_b.ptr(), Mrhs_b);
      rhs_B[ii] = rhs_b;
    }
  }

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

  // Build extended matrix
  this->phiLinSolv_->buildK(p);
  this->phiLinSolv_->buildb(rhs_B);
  Atilde_ = this->phiLinSolv_->buildATilde(cdt);

  // Build initial vector and compute matrix exponential in place
  auto v = this->phiLinSolv_->buildv(Atilde_->domain(), rhs_B[0]);
  Thyra::SolveStatus<Scalar> sStatus = this->matrixExponential(Atilde_, v);

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  //Atilde_->describe(*out, Teuchos::VERB_EXTREME);
  //Atilde_->domain()->describe(*out, Teuchos::VERB_EXTREME);

  // Get the first block of the multi-vector calculated from 2x2 multi-matrix
  auto v0 = v->getVectorBlock(0);  // V block
  Thyra::copy(*v0, x.ptr());

  return sStatus;
}

template <class Scalar>
Thyra::SolveStatus<Scalar>
PhiEvaluatorTaylor<Scalar>::matrixExponential(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
					      const Teuchos::RCP<Thyra::VectorBase<Scalar>> v)
{
  const int expansionOrder = getTaylorExpansionOrder();

  TEUCHOS_TEST_FOR_EXCEPTION(
      expansionOrder < 0,
      std::invalid_argument,
      "matrixExponential: expansionOrder must be nonnegative");

  // exp(L) * v is in range(L)
  const auto rangeSpace = L->range();

  // Iteration vector d_0 = v
  Teuchos::RCP<Thyra::VectorBase<Scalar>> d_k = Thyra::createMember(rangeSpace);
  Thyra::assign(d_k.ptr(), *v);

  // v = d_0 / 0! = d_0

  // allocate temporary vector
  Teuchos::RCP<Thyra::VectorBase<Scalar>> next = Thyra::createMember(rangeSpace);

  Thyra::SolveStatus<Scalar> sStatus;
  Scalar norm_d_k;
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
    Thyra::Vp_V(v.ptr(), *d_k);

    norm_d_k = Thyra::norm_inf(*d_k);

    // overflow is an upper bound on Thyra::norm_inf(*v);
    // it tracks how large the update matExpTemp could have gotten in intermediate iterations
    // any number larger than overflow * machine_eps should not be affected much by roundoff
    overflow += norm_d_k;

    // TODO: refine this and make dependent on Scalar type
    const Scalar cutoff = 1e22;
    if (overflow > cutoff)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
      break;
    }

    // terminate if the update drops below likely significance
    if (norm_d_k < overflow / cutoff)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }

    // set status if expansionOrder has been reached
    if (k == expansionOrder)
    {
      sStatus.achievedTol = norm_d_k;
      //sStatus.solveStatus = Thyra::SOLVE_STATUS_UNKNOWN;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }
  }

  std::stringstream ss;
  ss << "Taylor: Norm of solution=" << Thyra::norm_inf(*v)
     << " overflow=" << overflow
     << " final update=" << norm_d_k
     << " achieved in it. " << k << ".";
  sStatus.message = ss.str();

  std::cout << sStatus.message << std::endl;

  return sStatus;
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  setTaylorExpansionOrder(pl->get<int>("Taylor Expansion Order", 10));

  std::cout << "\nParameter List: " << *pl << std::endl;
  std::cout << "Taylor Expansion Order is " << getTaylorExpansionOrder() << std::endl;
}


// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> phi = Teuchos::rcp(new PhiEvaluatorTaylor<Scalar>());
  phi->setName("From createPhiEvaluatorTaylor");

  if (pl != Teuchos::null)
    phi->setPhiEvaluatorValues(pl);

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorTaylor_impl_hpp

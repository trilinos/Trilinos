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
Thyra::SolveStatus<Scalar> PhiEvaluatorTaylor<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
								  int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mrhs_b)
{
  // phi->setLumpMassMatrix(useLumpedMass_);
  this->phiLinSolv_->setLumpMassMatrix(this->lumpMassMatrix_);
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
    Thyra::Vp_V(matExpTemp.ptr(), *d_k);

    norm_d_k = Thyra::norm_inf(*d_k);

    // overflow is an upper bound on Thyra::norm_inf(*matExpTemp);
    // it tracks how large the update matExpTemp could have gotten in intermediate iterations
    // any number larger than overflow * machine_eps should not be affected much by roundoff
    overflow += norm_d_k;
      
    //std::cout << "Norm in it " << k << " of: solution=" << Thyra::norm_inf(*matExpTemp)
    //	      << " overflow=" << overflow << " update=" << norm_d_k << std::endl;

    // terminate if the update drops below likely significance
    if (norm_d_k < 1e-21 * overflow)
      break;
  }
  std::cout << "Taylor: Norm in it " << k << " of: solution=" << Thyra::norm_inf(*matExpTemp)
	    << " overflow=" << overflow << " update=" << norm_d_k << std::endl;

  //TODO return overflow * 1e-17 as an error estimate.

  matExp_v_ = matExpTemp; // This is required to wrap multivector as linearop

  return matExp_v_;
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

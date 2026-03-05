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
#include "Teuchos_Assert.hpp"
#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

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
      "Max Leja Order", 500,
      "Maximal order of the Leja polynomial used.\n"
      "\n"
      "The default is 500.");

  pl->set<double>(
      "leja_tol", 1.0e-18,
      "Leja polynomial convergence tolerance. Default is 1e-18");

  pl->set<double>(
      "leja_a", -1.0,
      "Minimum real bound of ellipse bounding system spectrum. The default is -1.0.");

  pl->set<double>(
      "leja_b", 0.0,
      "Maximum real bound of ellipse bounding system spectrum. The default is 0.0.");

  pl->set<double>(
      "leja_c", 0.5,
      "Maximum complex bound of the  ellipse bounding system spectrum. The default is 0.5.");

  return pl;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorLeja<Scalar>::computeLinOpPhi(const int phi_order,
                     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
                     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v)
{
  // TODO: why not use uint instead of int then?
  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order < 0,
      std::invalid_argument,
      "LinOpPhi: phi_order must be nonnegative.");

  const int expansionOrder = this->maxLejaOrder_;

  // TODO: optional fractional step size exp(tau*A_tilde)*v
  const Scalar tau = 1.0;

  // phi_k(L) * v is in range(L)
  const auto rangeSpace = L->range();

  // compute shift and scale parameters
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();

  // TODO: update the divided differences (or read from cache)

  // Iteration vector vm_0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> vm_k = Thyra::createMember(rangeSpace);
  // Iteration vector qm_0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> qm_k = Thyra::createMember(rangeSpace);
  // Temp storage for Matvec result
  Teuchos::RCP<Thyra::VectorBase<Scalar>> av = Thyra::createMember(rangeSpace);

  // 0th term of the leja polynomial
  std::complex<double> coeff = this->lp_dd_.at(0);
  Scalar coeff_re = Scalar(coeff.real());
  Thyra::assign(vm_k.ptr(), *v);  // w_m
  Thyra::V_StV(v, coeff_re, *vm_k);

  // storage for error est
  Scalar norm_vm_k = Thyra::norm_inf(*vm_k);
  Scalar overflow = 0.;
  const Scalar cutoff = 1e22;
  Thyra::SolveStatus<Scalar> sStatus;

  Scalar lp_sc_re;
  Scalar lp_sc_im;

  int k = 1;
  while (k < expansionOrder-1)
  {
    // extract divided diff
    std::complex<double> coeff = this->lp_dd_.at(k);
    coeff_re = Scalar(coeff.real());

    // Real leja point case
    if (this->lp_base_.at(k-1).lt == LejaType::LPREAL) {
      // compute shifted and scaled leja point
      lp_sc_re = Scalar( shift + scale * this->lp_base_.at(k-1).get().at(0).real() );
      // av = (tau*A)*vm
      Thyra::apply(*L, Thyra::NOTRANS, *vm_k, av.ptr(), tau, 1.0);
      // vm_k = (av - lp_re[k-1]*vm_k)
      Thyra::V_VpStV(vm_k.ptr(), *av, -lp_sc_re, *vm_k);
      // vm_k = vm_k / scale
      Thyra::V_StV(vm_k.ptr(), 1.0 / scale, *vm_k);
      // add vm_k*coeff to the final result
      Thyra::Vp_StV(v, coeff_re, *vm_k);
      k += 1;

      norm_vm_k = Thyra::norm_inf(*vm_k) * coeff_re;
    }
    else if (this->lp_base_.at(k-1).lt == LejaType::LPCONJ)  {
      // first update
      lp_sc_re = Scalar( shift + scale * this->lp_base_.at(k-1).get().at(0).real() );
      Thyra::apply(*L, Thyra::NOTRANS, *vm_k, av.ptr(), tau, 1.0);
      Thyra::V_VpStV(qm_k.ptr(), *av, -lp_sc_re, *vm_k);
      Thyra::V_StV(qm_k.ptr(), 1.0 / scale, *qm_k);
      Thyra::Vp_StV(v, coeff_re, *qm_k);
      k += 1;

      // conjugate update
      lp_sc_re = Scalar( shift + scale * this->lp_base_.at(k-1).get().at(1).real() );
      lp_sc_im = Scalar( shift + scale * this->lp_base_.at(k-1).get().at(1).imag() );
      std::complex<double> coeff = this->lp_dd_.at(k+1);
      coeff_re = Scalar(coeff.real());
      Thyra::apply(*L, Thyra::NOTRANS, *qm_k, av.ptr(), tau, 1.0);
      Thyra::V_VpStV(vm_k.ptr(), *av, -lp_sc_re, *qm_k);
      Thyra::V_StV(vm_k.ptr(), 1.0 / scale, *vm_k);
      Thyra::Vp_StV(vm_k.ptr(), (lp_sc_im / scale) * (lp_sc_im / scale), *vm_k);
      Thyra::Vp_StV(v, coeff_re, *vm_k);
      k += 1;

      norm_vm_k = Thyra::norm_inf(*vm_k) * coeff_re;
    }
    else {
      // TODO: ERROR
      TEUCHOS_ASSERT(false);
    }
    overflow += norm_vm_k;

    // TODO: refine this and make dependent on Scalar type
    const Scalar cutoff = 1e22;
    if (overflow > cutoff)
    {
      sStatus.achievedTol = norm_vm_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
      break;
    }

    // terminate if the update drops below likely significance
    if (norm_vm_k < overflow / cutoff)
    {
      sStatus.achievedTol = norm_vm_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }

    // terminate if the update drops below user tol
    if (norm_vm_k < this->leja_tol_)
    {
      sStatus.achievedTol = norm_vm_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }
  }

  std::stringstream ss;
  ss << "Leja: Norm of solution=" << Thyra::norm_inf(*v)
     << " overflow=" << overflow
     << " final update=" << norm_vm_k
     << " achieved in it. " << k << ".";
  sStatus.message = ss.str();

  return sStatus;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  maxLejaOrder_ = pl->get<int>("Max Leja Order", 500);
  leja_tol_ = pl->get<double>("leja_tol", 1.0e-18);
  leja_a_ = pl->get<double>("leja_a", -1.0);
  leja_b_ = pl->get<double>("leja_b", 0.0);
  leja_c_ = pl->get<double>("leja_c", 0.5);

  // TODO: has to be set to true, only matrix exponential is implemented
  this->useAtildeForSingleRHS_ = true;

  std::cout << "\nuseAtildeForSingleRHS_: " << this->useAtildeForSingleRHS_ << std::endl;
  std::cout << "Parameter List: " << *pl << std::endl;
  std::cout << "Leja Order is " << maxLejaOrder_ << std::endl;
}

template <class Scalar>
std::tuple<double, double> PhiEvaluatorLeja<Scalar>::getShiftScale()
{
  // real half axis
  double hx_re = (leja_b_ - leja_a_) / 2.0;
  // imaj half axis
  double hx_im = leja_c_;
  // leja ellipse shift and scale parameters
  double shift = (leja_a_ + leja_b_) / 2.0;
  double scale = (hx_re + hx_im) / 2.0;
  return std::make_tuple(shift, scale);
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

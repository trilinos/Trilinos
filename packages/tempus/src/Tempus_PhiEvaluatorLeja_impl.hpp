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

#include <cmath>
#include <complex>
#include "Tempus_PhiEvaluatorLeja.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluatorLeja_decl.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_ArrayRCPDecl.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCPNode.hpp"
#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorStdOps_decl.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "TpetraCore_config.h"
#include "Tpetra_CombineMode.hpp"

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
      "Expansion Order", 300,
      "Order of the Leja polynomial used.\n"
      "\n"
      "The default is 300.");

  pl->set<int>(
      "Leja DD Method", 1,
      "DD Method to use. 0 for Recurrence. 1 for Taylor Series. 2 for DD_phi. 3 for Taylor Real.");

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

  pl->set<double>(
      "Leja Ellipse Saftey Factor", 1.0,
      "Saftey scaling factor applied to the adapted Leja ellipse size. The default is 1.0.");

  return pl;
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorLeja<Scalar>::computeLinOpPhi(const int phi_order,
                     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
                     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
                     const Scalar cdt)
{
  // use ScalarTraits instead of std::abs (magnitude, returns a magnitudeType)
  using ST = typename Teuchos::ScalarTraits<Scalar>;
  using magScalar = typename ST::magnitudeType;

  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order != 0,
      std::invalid_argument,
      "LinOpPhi: phi_order must be zero.");

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor phitimer(*timerPhi_);
#endif

  // phi_k(L) * v is in range(L)
  const auto rangeSpace = L->range();

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  // scale the Leja ellipse by the provided (fractional) timestep
  const Scalar scale = cdt * Scalar(std::get<0>(transform_params));

  // Get lejaOrder divided differences, we need one more than the expansion order (polynomial order)
  const int lejaOrder = this->getExpansionOrder() + 1;
  auto lp_dd = getDividedDiffs(phi_order, cdt, lejaOrder);
  TEUCHOS_ASSERT(lp_dd.size() == lejaOrder);

  //std::cout << "DD: " << std::endl;
  //for (const auto& dd : lp_dd) {
  //  std::cout << dd << ' ';
  //}
  //std::cout << std::endl;

  // Iteration vector vm_0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> vm_k = Thyra::createMember(rangeSpace);
  // Iteration vector qm_0
  Teuchos::RCP<Thyra::VectorBase<Scalar>> qm_k = Thyra::createMember(rangeSpace);
  // Temp storage for Matvec result
  Teuchos::RCP<Thyra::VectorBase<Scalar>> av = Thyra::createMember(rangeSpace);

  // 0th term of the leja polynomial
  Scalar coeff = lp_dd[0];
  //std::cout << "c[0]: " << coeff << std::endl;

  Thyra::V_V(vm_k.ptr(), *v);
  Thyra::V_StV(v, coeff, *vm_k);

  // storage for error est
  magScalar norm_vm_k = Thyra::norm_inf(*vm_k);
  // norm of the update
  magScalar norm_d_k = ST::magnitude(coeff) * norm_vm_k;
  // upper bound on solution size
  magScalar overflow = norm_d_k;
  Thyra::SolveStatus<Scalar> sStatus;

  // leja point index k starts at 1
  int k = 1;
  // leja polynomial term index lp_k starts at 0
  for (int lp_k = 0; k < lejaOrder && lp_k < lejaPointsBase_.size(); k++, lp_k++)
  {
    // print the update vector vm_k
    //v->describe(*this->getOStream(), Teuchos::VERB_EXTREME);
    //std::cout << "Norm d_k: " << norm_d_k << " v_k: " << norm_vm_k << std::endl;

    // compute transformed unscaled Leja point
    const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_k], transform_params);

    const Scalar lp_re = Scalar( lp.lp.real() );
    const Scalar lp_im = Scalar( lp.lp.imag() );

    // Real leja point case
    if (lp.lpt == LpType::LPREAL) {
      // extract divided diff
      coeff = lp_dd[k];
      // std::cout << "c,lp,shift,scale: " << coeff << " " << lp_sc.lp << " real " << shift << " " << scale << std::endl;

      {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor linoptimer(*timerLinOp_);
#endif
        // copy vm_k to temp vector
        Thyra::V_V(av.ptr(), *vm_k);
        // vm_k = (L@vm_k - lp_re[k-1]*vm_k) / scale
        Thyra::apply(*L, Thyra::NOTRANS, *av, vm_k.ptr(), 1 / scale, -lp_re);
      }

      // add vm_k*coeff to the final result
      Thyra::Vp_StV(v, coeff, *vm_k);

      norm_vm_k = Thyra::norm_inf(*vm_k);
      norm_d_k = ST::magnitude(coeff) * norm_vm_k;
      overflow += norm_d_k;
    }
    else if (lp.lpt == LpType::LPCONJ) {
      // extract divided diff
      coeff = lp_dd[k];
      //std::cout << "c,lp,shift,scale: " << coeff << " " << lp_sc.lp << " " << shift << " " << scale << std::endl;

      // copy vm_k to qm_k vector to save it
      Thyra::V_V(qm_k.ptr(), *vm_k);

      {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor linoptimer(*timerLinOp_);
#endif
        // vm_k = (L@vm_k - lp_re*vm_k) / scale
        Thyra::apply(*L, Thyra::NOTRANS, *qm_k, vm_k.ptr(), 1 / scale, -lp_re);
      }

      // add vm_k*coeff to the final result
      Thyra::Vp_StV(v, coeff, *vm_k);
      norm_vm_k = Thyra::norm_inf(*vm_k);
      norm_d_k = ST::magnitude(coeff) * norm_vm_k;
      overflow += norm_d_k;

      // increment polynomial degree, but keep Leja point and handle conjugate pair
      k++;
      // std::cout << "qm" << std::endl;
      // vm_k->describe(*this->getOStream(), Teuchos::VERB_EXTREME);

      if (k < lp_dd.size())
      {
        // conjugate update
        coeff = lp_dd[k];

        //std::cout << "Norm d_k: " << norm_d_k << " v_k: " << norm_vm_k << std::endl;
        //std::cout << "c,lp: " << coeff << " " << std::conj(lp_sc.lp) << std::endl;

        {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor linoptimer(*timerLinOp_);
#endif
          // copy vm_k to a new temp vector (don't overwrite qm_k)
          Thyra::V_V(av.ptr(), *vm_k);
          // vm_k = (L@vm_k - lp_re*vm_k) / scale + ((lp_im/scale)**2)*qm_k
          Thyra::apply(*L, Thyra::NOTRANS, *av, vm_k.ptr(), 1/scale, -lp_re);
          Thyra::Vp_StV(vm_k.ptr(), lp_im * lp_im, *qm_k);
        }

        // add vm_k*coeff to the final result
        Thyra::Vp_StV(v, coeff, *vm_k);
        norm_vm_k = Thyra::norm_inf(*vm_k);
        norm_d_k = ST::magnitude(coeff) * norm_vm_k;
        overflow += norm_d_k;
      }
    }
    else {
      // We do not support non-conjugate complex Leja points
      TEUCHOS_ASSERT(false);
    }

    // TODO: refine this and make dependent on Scalar type
    const Scalar cutoff = 1e22;
    if (overflow > cutoff)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_UNCONVERGED;
      break;
    }

    // terminate if the update drops below user tol
    if (k >= lejaOrder || norm_d_k < this->leja_tol_)
    {
      sStatus.achievedTol = norm_d_k;
      sStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
      break;
    }
  }

  std::stringstream ss;
  ss << "Leja: Norm of solution=" << Thyra::norm_inf(*v)
     << " overflow=" << overflow
     << " final update=" << norm_d_k
     << " iteration vector=" << norm_vm_k
     << " achieved in it. " << k << ".";
  sStatus.message = ss.str();

  // std::cout << sStatus.message << std::endl;

  return sStatus;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::initLejaPointsBase(const int lejaOrder)
{
  // always have two points at least
  int maxLejaOrder = std::max(2, lejaOrder);

  if (lejaOrder < 0)
  {
    // The maximum number of Leja points is: lejaOrder + 1, but due to conjugacy, we technically need less points stored:
    // For real Leja points given as the real part of the points below,
    //    this is perfect, since every conjugate Leja point maps to a single real point
    // For imaginary Leja points given as the imaginary part,
    //    this mostly correct but sufficient for lejaOrder + 1 > 2, since the first two points map both to 0
    // For conjugate Leja points, we only need 2 + lejaOrder / 2,
    //    this still works since it is an upper bound
    maxLejaOrder = std::max(2, this->getExpansionOrder() + 1);
  }

  lejaPointsBase_ = Teuchos::arcp<LejaPoint>(maxLejaOrder);

  lejaPointsBase_[0] = {-1., LPREAL};
  lejaPointsBase_[1] = {1., LPREAL};

  std::complex<double> root_unity(0, 1);
  int full_half_circle = 1;
  for (int lpk = 2; lpk < lejaPointsBase_.size(); lpk++)
  {
    // get the old leja Point from the last full half circle
    std::complex<double> next_lp = lejaPointsBase_[lpk - full_half_circle].lp;

    // rotate the leja point by root of unity
    next_lp *= root_unity;

    // save the new leja point
    lejaPointsBase_[lpk] = {next_lp, LPCONJ};

    // if we have completed one full half circle (upper complex half plane)
    if (lpk >= 2*full_half_circle)
    {
      full_half_circle *= 2;
      root_unity = std::sqrt(root_unity);
    }
  }

  // swap the first two real leja points (to have 1 first, not essential)
  std::swap(lejaPointsBase_[0], lejaPointsBase_[1]);
}

template <class Scalar>
constexpr std::tuple<double, double, double> PhiEvaluatorLeja<Scalar>::getScaleFromBase()
{
  // real half axis
  const double scale_re = (leja_b_ - leja_a_) / 2.0;
  // imag half axis
  const double scale_im = leja_c_;
  // Leja ellipse shift and scale parameters
  const double scale = (scale_im + scale_re) / 2.0;
  const double shift = (leja_a_ + leja_b_) / 2.0;

  // normalize the transform parameters for later convenience
  const double sc_shift = shift/scale;
  const double sc_anis = scale_re/scale; // in [0, 2] interval

  return std::make_tuple(scale, sc_shift, sc_anis);
}

template <class Scalar>
constexpr LejaPoint PhiEvaluatorLeja<Scalar>::transformLejaPoint(const LejaPoint& lp_base,
                                                                 const std::tuple<const double, const double, const double>& scale_params)
{
  // copy base Leja point
  LejaPoint scaled_lp = lp_base;

  const double sc_shift = std::get<1>(scale_params), sc_anis = std::get<2>(scale_params);

  // apply normalized transform, not scaling
  scaled_lp.lp.real(sc_anis * scaled_lp.lp.real() + sc_shift);
  scaled_lp.lp.imag((2.0 - sc_anis) * scaled_lp.lp.imag());

  return scaled_lp;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setExpansionOrder(const int expansionOrder)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      expansionOrder <= 0,
      std::invalid_argument,
      "setExpansionOrder: order must be positive.");

  expansionOrder_ = expansionOrder;

  initLejaPointsBase();
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setLejaEllipse(const double a, const double b, const double c)
{
  TEUCHOS_ASSERT(a <= b);
  TEUCHOS_ASSERT(c >= 0.0);
  leja_a_ = a;
  leja_b_ = b;
  leja_c_ = c;
}

template <class Scalar>
Teuchos::Tuple<double, 3> PhiEvaluatorLeja<Scalar>::getLejaEllipse()
{
  return Teuchos::tuple<double>(leja_a_, leja_b_, leja_c_);
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::adaptEvaluator()
{
  this->phiLinSolv_->computeJacobianSpectrumBounds(leja_a_, leja_b_, leja_c_);
  // scale ellipse by saftey factor
  leja_a_ *= leja_sf_;
  leja_c_ *= leja_sf_;
  std::cout << "Adapted Leja ellipse parameters. a=" <<
    leja_a_ << " b=" << leja_b_ << " c=" << leja_c_ << std::endl;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setDivideDifferenceMethod(const int ddMethod)
{
  ddMethod_ = ddMethod;
}

template <class Scalar>
LejaPoint PhiEvaluatorLeja<Scalar>::getLpSc(const int lp_idx)
{
  TEUCHOS_ASSERT(lp_idx < lejaPointsBase_.size());

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  const double scale = std::get<0>(transform_params);

  // transform but do not scale Leja point
  const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_idx], transform_params);

  // return scaled and transformed Leja point
  LejaPoint lp_sc = LejaPoint{scale * lp.lp, lp.lpt};
  return lp_sc;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  leja_sf_ = pl->get<double>("Leja Ellipse Saftey Factor", 1.0);
  leja_tol_ = pl->get<double>("leja_tol", 1.0e-18);
  setDivideDifferenceMethod(pl->get<int>("Leja DD Method", 2));
  setExpansionOrder(pl->get<int>("Expansion Order", 300));

  // TODO: has to be set to true, only matrix exponential is implemented
  this->useAtildeForSingleRHS_ = true;

  std::cout << "\nuseAtildeForSingleRHS_: " << this->useAtildeForSingleRHS_ << std::endl;
  std::cout << "Parameter List: " << *pl << std::endl;
  std::cout << "Expansion Order is " << getExpansionOrder() << std::endl;

  setLejaEllipse(
    pl->get<double>("leja_a", -1.0),
    pl->get<double>("leja_b", 0.0),
    pl->get<double>("leja_c", 0.5)
  );
}

template <class Scalar>
Teuchos::ArrayRCP<Scalar> PhiEvaluatorLeja<Scalar>::getDividedDiffs(const int k, const Scalar cdt, const int lejaOrder)
{
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor ddtimer(*timerDD_);
#endif

  // ensure that we have enough Leja points
  if (lejaOrder > lejaPointsBase_.size()) {
    initLejaPointsBase(lejaOrder);
  }

  switch (ddMethod_)
  {
  case 0:
    return getDividedDiffsRC(k, cdt, lejaOrder);
    break;
  case 1:
    return getDividedDiffsTS(k, cdt, lejaOrder);
    break;
  case 2:
    return getDividedDiffsPhi(k, cdt, lejaOrder);
    break;
  case 3:
  default:
    return getDividedDiffsTSR(k, cdt, lejaOrder);
  }
}

template <class Scalar>
Teuchos::ArrayRCP<Scalar> PhiEvaluatorLeja<Scalar>::getDividedDiffsRC(const int phi_order, const Scalar cdt, const int lejaOrder)
{
  TEUCHOS_ASSERT(phi_order == 0);

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  // scale the Leja ellipse by the provided (fractional) timestep
  const double scale = (double)cdt * std::get<0>(transform_params);
  // since leja points use double, shift should too for this method

  Teuchos::ArrayRCP<std::complex<double>> x = Teuchos::arcp<std::complex<double>>(lejaOrder);
  Teuchos::ArrayRCP<std::complex<double>> d_x = Teuchos::arcp<std::complex<double>>(lejaOrder);

  // initialize list of Leja points and function values
  for (int idx = 0, lp_idx = 0; idx < lejaOrder && lp_idx < lejaPointsBase_.size(); idx++, lp_idx++)
  {
    const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_idx], transform_params);
    if (lp.lpt == LPCONJ)
    {
      x[idx] = lp.lp;
      d_x[idx] = std::exp(scale * x[idx]);
      if (++idx < lejaOrder)
      {
        x[idx] = std::conj(lp.lp);
        d_x[idx] = std::exp(scale * x[idx]);
      }
    }
    else
    {
      x[idx] = lp.lp;
      d_x[idx] = std::exp(scale * x[idx]);
    }
  }

  for (int idx = 0; idx < lejaOrder-1; idx++)
  {
    // Compute the next set of divided differences
    for (int idy = idx+1; idy < lejaOrder; idy++)
    {
      d_x[idy] = (d_x[idy] - d_x[idx]) / (x[idy] - x[idx]);
    }
  }

  // read out real part and convert to desired Scalar type
  Teuchos::ArrayRCP<Scalar> dd_phi = Teuchos::arcp<Scalar>(lejaOrder);
  for (int idx = 0; idx < lejaOrder - 1; idx++)
    dd_phi[idx] = Scalar(d_x[idx].real());

  return dd_phi;
}

template <class Scalar>
Teuchos::ArrayRCP<Scalar> PhiEvaluatorLeja<Scalar>::getDividedDiffsPhi(
    const int phi_order, const Scalar cdt, const int lejaOrder)
{
  // Compute Newton divided differences of phi_{phi_order} at the Leja points using
  // the Zivcovich (2019) H-factorization + scaling-and-squaring method.
  //
  // Ref: F. Zivcovich. "Fast and accurate computation of divided differences for
  //      analytic functions, with an application to the exponential function."
  //      Dolomites Research Notes on Approximation. 12. 2019.

#ifdef HAVE_TEUCHOS_COMPLEX
  using cplx = std::complex<double>;

  const int n_leja  = lejaOrder;                   // number of output coefficients
  const int n_total = phi_order + n_leja;          // total interpolation points
  const int p       = 30;                          // Taylor truncation terms
  const int cap_n   = n_total - 1 + p;             // Taylor degree N in the paper

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  // scale the Leja ellipse by the provided (fractional) timestep
  const double scale = (double)cdt * std::get<0>(transform_params);

  // build z[0..n_total]
  //  in contrast to Zivcovich we only shift the Leja points, but not scale them.
  //  this avoids over and underflow, since the original leja points are of roghly unit size
  // z[0..phi_order]         = 0  (leading zeros encode the phi_k recurrence)
  // z[phi_order..n_total]     = shift/scale + lp[i], expanding LPCONJ pairs to (lp, conj(lp))
  Teuchos::ArrayRCP<cplx> z = Teuchos::arcp<cplx>(n_total);
  {
    int z_idx = 0;
    for ( ; z_idx < phi_order; z_idx++)
      z[z_idx] = cplx(0.0);
    for (int lp_idx = 0; lp_idx < lejaPointsBase_.size() && z_idx < n_total; ++lp_idx)
    {
      const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_idx], transform_params);
      if (lp.lpt == LPCONJ)
      {
        // upper-half-plane point
        z[z_idx++] = lp.lp;
        // conjugate (lower-half-plane)
        if (z_idx < n_total)
          z[z_idx++] = std::conj(lp.lp);
      }
      else  // LPREAL: zero out any floating-point imaginary noise
      {
        z[z_idx++] = cplx(lp.lp.real(), 0.0);
      }
    }
    // ensure we did not run out of points in lejaPointsBase_
    TEUCHOS_ASSERT(z_idx == n_total);
  }

  // for optimal alignment around center, shift again by mean mu, will be corrected at the end
  cplx mu(0.0, 0.0);
  for (int i = 0; i < n_total; ++i)
    mu += z[i];
  mu /= double(n_total);
  for (int i = 0; i < n_total; ++i)
    z[i] -= mu;

  // build lower-triangle of F, containing the pairwise differences
  //  this is just a performance optimization (does it even help?)
  //  in contrast to Zivcovich, scaled points are taken, see above
  Teuchos::SerialDenseMatrix<int, cplx> F_mat(n_total, n_total);
  F_mat.putScalar(cplx(0.0, 0.0));
  for (int i0 = 0; i0 < n_total - 1; ++i0)
    for (int j0 = i0 + 1; j0 < n_total; ++j0)
      F_mat(j0, i0) = z[i0] - z[j0];

  // number of squarings (maximum distance of any given point)
  double max_abs = 0.0;
  for (int i0 = 0; i0 < n_total - 1; ++i0)
    for (int j0 = i0 + 1; j0 < n_total; ++j0)
    {
      double v = std::abs(F_mat(j0, i0));
      if (v > max_abs) max_abs = v;
    }

  // correct for the scaling difference to Zivkovich paper
  mu *= scale;
  max_abs *= scale;

  // Find the number of subdivisions s of the unit interval such that max(F) / s < 3.5
  //   the corresponding number of squarings is log2(s), if s happens to be a power of 2
  // alternative:
  // const int log2_s   = std::max(int( std::ceil(std::log2(max_abs / 3.5)), 0);
  // const double s_dbl = spd::pow(double(log2_s), 2.0);
  const int    s     = std::max(int(std::ceil(max_abs / 3.5)), 1);
  const double s_dbl = double(s);

  // seed dd[0..cap_n]: dd[kk] = scale^kk / (kk! * s^kk)
  //  the factor scale^kk is added to correct for the scaling difference (this mitigates underflow if scale and s are large)
  Teuchos::ArrayRCP<cplx> dd = Teuchos::arcp<cplx>(cap_n + 1);
  std::fill(dd.begin(), dd.end(), cplx(0.0, 0.0));
  dd[0] = cplx(1.0, 0.0);
  // avoid overflow but tolerate underflow
  double running_fraction = 1.0;
  // std::cout << "initial dd:  ";
  for (int kk = 1; kk <= cap_n; ++kk)
  {
    running_fraction *= scale / (kk * s_dbl);
    dd[kk] = cplx(running_fraction, 0.0);
    // std::cout << dd[kk] << ", ";
  }
  // std::cout << std::endl;

  // H-factorization sweep:
  //  In contrast to Zivcovich, the z points are scaled.
  //  Thus our dd[j] and the columns of F_mat (the upper triangle part, not the lower z_i-z_j precomputed points)
  //  contain an additional factor scale^k, which we want.
  for (int j = n_total - 1; j >= 0; --j)
  {
    // First inner loop: Taylor remainder sweep
    // k0 from cap_n-1 down to n_total-2-j (exclusive): dd[k0] += z[j] * dd[k0+1]
    for (int k0 = cap_n - 1; k0 > n_total - 2 - j; --k0)
      dd[k0] = dd[k0] + z[j] * dd[k0 + 1];

    // Second inner loop: divided-difference sweep using F lower triangle
    // k0 from (n_total-2-j) downto 0: dd[k0] += F(k0+j+1, j) * dd[k0+1]
    for (int k0 = n_total - 2 - j; k0 >= 0; --k0)
      dd[k0] = dd[k0] + F_mat(k0 + j + 1, j) * dd[k0 + 1];

    // Store dd[0..n_total-j] into upper-triangle row j of F
    for (int col = 0; col < n_total - j; ++col)
      F_mat(j, j + col) = dd[col];
  }

  // overwrite diagonal: F(i,i) = exp(scale * z[i] / s)
  for (int i = 0; i < n_total; ++i)
    F_mat(i, i) = std::exp(cplx(scale / s_dbl, 0.0) * z[i]);

  // zero lower triangle to remove the precomputed z_i-z_j
  for (int i = 1; i < n_total; ++i)
    for (int j = 0; j < i; ++j)
      F_mat(i, j) = cplx(0.0, 0.0);

  // substepping s times (squaring)
  // build the vector that is F_mat * unit_vector_e1
  Teuchos::SerialDenseMatrix<int, cplx> dd_row(1, n_total);
  for (int j = 0; j < n_total; ++j)
    dd_row(0, j) = F_mat(0, j);

  // apply the matrix s-1 times to the initial vector.
  Teuchos::SerialDenseMatrix<int, cplx> tmp_row(1, n_total);
  for (int i_s = 0; i_s < s - 1; ++i_s)
  {
    // copy dd_row for matvec
    tmp_row = dd_row;
    // dd_row = 1.0 * tmp_row * F_mat + 0.0 * dd_row
    dd_row.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, cplx(1.0, 0.0), tmp_row, F_mat, cplx(0.0, 0.0));
  }

  // out[i] = exp(mu) * dd_row[phi_order + i]
  // use only the real part and convert to Scalar
  Teuchos::ArrayRCP<Scalar> dd_phi = Teuchos::arcp<Scalar>(n_leja);
  const cplx exp_mu = std::exp(mu);
  // std::cout << "final dd:  ";
  for (int i = 0; i < n_leja; ++i)
  {
    const cplx dd_i = exp_mu * dd_row(0, phi_order + i);
    dd_phi[i] = Scalar(dd_i.real());
    // std::cout << dd_i << ", ";
  }
  // std::cout << std::endl;

  return dd_phi;
#else
  std::cout << "WARNING: getDividedDiffsPhi requires Trilinos configured with Trilinos_ENABLE_COMPLEX=ON" << std::endl;
  std::cout << "WARNING: falling back to getDividedDiffsTSR." << std::endl;
  return getDividedDiffsTSR(phi_order, cdt, lejaOrder);
#endif
}

template <class Scalar>
Teuchos::ArrayRCP<Scalar> PhiEvaluatorLeja<Scalar>::getDividedDiffsTS(const int phi_order, const Scalar cdt, const int lejaOrder)
{
  TEUCHOS_ASSERT(phi_order == 0);

#ifdef HAVE_TEUCHOS_COMPLEX
  using cplx = std::complex<double>;

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  // scale the Leja ellipse by the provided (fractional) timestep
  const double scale = (double)cdt * std::get<0>(transform_params);

  // build the shifted and scaled Hm matrix
  Teuchos::SerialDenseMatrix<int, cplx> Hm(lejaOrder, lejaOrder);
  // diagonal elements are the leja points
  int dd_idx = 0;
  int lp_idx = 0;
  while (dd_idx < lejaOrder) {
    const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_idx], transform_params);

    // conj lp case
    if (lp.lpt == LPCONJ) {
      if (dd_idx == lejaOrder) break;
      Hm(dd_idx, dd_idx) = scale * lp.lp;
      if (dd_idx+1 < lejaOrder) Hm(dd_idx+1, dd_idx) = scale;
      dd_idx += 1;
      if (dd_idx == lejaOrder) break;
      Hm(dd_idx, dd_idx) = scale * std::conj(lp.lp);
      if (dd_idx+1 < lejaOrder) Hm(dd_idx+1, dd_idx) = scale;
      dd_idx += 1;
    }
    else {
      if (dd_idx == lejaOrder) break;
      Hm(dd_idx, dd_idx) = scale * lp.lp;
      if (dd_idx+1 < lejaOrder) Hm(dd_idx+1, dd_idx) = scale;
      dd_idx += 1;
    }
    lp_idx += 1;
  }

  //std::cout << "LP: " << std::endl;
  //for (const auto& lp : lp_) {
  //  std::cout << lp.lp << ' ';
  //}
  //std::cout << std::endl;

  //std::cout << "Hm: " << std::endl;
  //Hm.print(std::cout);

  // compute diagonal mean
  cplx diag_sum = cplx(0.0, 0.0);
  for (int i=0; i < lejaOrder; ++i) {
    diag_sum += Hm(i, i);
  }
  cplx mu = diag_sum / double(lejaOrder);

  // shift diagonal to zero mean
  for (int i=0; i < lejaOrder; ++i) {
    Hm(i, i) -= mu;
  }

  // Scaling
  double s_scale = Hm.normInf();
  int n_sq = std::max(int( std::ceil(std::log2(s_scale) - 1.) ), 1);

  double h_scale = 1.0 / std::pow(2.0, n_sq);
  Hm.scale(h_scale);

  // compute phi_0(Hm) by Taylor series
  //copy Hm to A
  Teuchos::SerialDenseMatrix<int, cplx> A(Teuchos::Copy, Hm);

  Teuchos::SerialDenseMatrix<int, cplx> Ts(lejaOrder, lejaOrder);
  Ts = 0.;

  for (int i=0; i < lejaOrder; ++i) {
    Ts(i, i) = std::complex(1.0, 0.0);
  }

  // Ts = I/(p!) + Hm^1/(1+p)! + Hm^2/(2+p)! ...
  int ts_order = 17;

  Teuchos::SerialDenseMatrix<int, cplx> Mtmp(lejaOrder, lejaOrder);
  Mtmp = 0.;

  for (int k = 1; k < ts_order; ++k) {
    Ts += A;
    // A = Hm^k/(k)!

    // Compute next A = Hm^(k+1)/(k+1)!
    double scale = 1. / (k+1);
    Mtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, scale, Hm, A, 0.0);
    A = Mtmp;
  }
  Ts += A;

  // Squaring
  for (int s=0; s < n_sq; ++s) {
    // TODO: Can this work without tmp output storage?
    Mtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ts, Ts, 0.0);
    Ts = Mtmp;
  }

  // unshift and extract first column
  // use only the real part and convert to Scalar
  Teuchos::ArrayRCP<Scalar> dd_phi = Teuchos::arcp<Scalar>(lejaOrder);
  const cplx exp_mu = std::exp(mu);
  for (int i = 0; i < lejaOrder; ++i)
  {
    const cplx dd_i = exp_mu * Ts(i, 0);
    dd_phi[i] = Scalar(dd_i.real());
  }
  return dd_phi;
#else
  std::cout << "WARNING: getDividedDiffsTS requires Trilinos configured with Trilinos_ENABLE_COMPLEX=ON" << std::endl;
  std::cout << "WARNING: falling back to getDividedDiffsTSR." << std::endl;
  return getDividedDiffsTSR(phi_order, cdt, lejaOrder);
#endif
}

template <class Scalar>
Teuchos::ArrayRCP<Scalar> PhiEvaluatorLeja<Scalar>::getDividedDiffsTSR(const int phi_order, const Scalar cdt, const int lejaOrder)
{
  TEUCHOS_ASSERT(phi_order == 0);

  // get scale and transform parameters
  std::tuple<double, double, double> transform_params = getScaleFromBase();
  // scale the Leja ellipse by the provided (fractional) timestep
  const Scalar scale = cdt * Scalar(std::get<0>(transform_params));

  // build the shifted and scaled Hm matrix
  Teuchos::SerialDenseMatrix<int, Scalar> Hm(lejaOrder, lejaOrder);
  // diagonal elements are the leja points

  for (int lp_idx = 0, dd_idx = 0; lp_idx < lejaPointsBase_.size() && dd_idx < lejaOrder; lp_idx++, dd_idx++) {
    const LejaPoint lp = transformLejaPoint(lejaPointsBase_[lp_idx], transform_params);
    Scalar lp_real = Scalar(lp.lp.real()), lp_imag = Scalar(lp.lp.imag());
    // conj lp case
    if (lp.lpt == LPCONJ) {
      Hm(dd_idx, dd_idx) = scale * lp_real;
      if (dd_idx + 1 < lejaOrder) Hm(dd_idx + 1, dd_idx) = scale;

      if (++dd_idx < lejaOrder) {
        Hm(dd_idx - 1, dd_idx) = - scale * lp_imag * lp_imag;
        Hm(dd_idx, dd_idx) = scale * lp_real;
        if (dd_idx + 1 < lejaOrder)
          Hm(dd_idx + 1, dd_idx) = scale;
      }
    }
    else {
      Hm(dd_idx, dd_idx) = scale * lp_real;
      if (dd_idx + 1 < lejaOrder)
        Hm(dd_idx + 1, dd_idx) = scale;
    }
  }

  // compute diagonal mean
  Scalar diag_sum = 0;
  for (int i = 0; i < lejaOrder; ++i) {
    diag_sum += Hm(i, i);
  }
  Scalar mu = diag_sum / Scalar(lejaOrder);

  // shift diagonal to zero mean
  for (int i = 0; i < lejaOrder; ++i) {
    Hm(i, i) -= mu;
  }

  // Scaling
  double s_scale = Hm.normInf();
  int n_sq       = std::max(int(std::ceil(std::log2(s_scale) - 1.)), 1);

  double h_scale = 1.0 / std::pow(2.0, n_sq);
  Hm.scale(Scalar(h_scale));

  // compute phi_0(Hm) by Taylor series
  //copy Hm to A
  Teuchos::SerialDenseMatrix<int, Scalar> A(Teuchos::Copy, Hm);

  Teuchos::SerialDenseMatrix<int, Scalar> Ts(lejaOrder, lejaOrder);
  Ts = 0.;

  for (int i = 0; i < lejaOrder; ++i) {
    Ts(i, i) = 1.0;
  }

  // Ts = I/(p!) + Hm^1/(1+p)! + Hm^2/(2+p)! ...
  int ts_order = 17;

  //std::cout << "ts_order: " << ts_order << std::endl;
  //std::cout << "s_scale: " << s_scale << std::endl;
  //std::cout << "mu: " << mu << std::endl;
  //std::cout << "n_sq: " << n_sq << std::endl;
  //std::cout << "h_scale: " << h_scale << std::endl;

  Teuchos::SerialDenseMatrix<int, Scalar> Mtmp(lejaOrder, lejaOrder);
  Mtmp = 0.;

  for (int k = 1; k < ts_order; ++k) {
    Ts += A;
    // A = Hm^k/(k)!

    // Compute next A = Hm^(k+1)/(k+1)!
    Scalar scale = Scalar(1. / (k+1));
    Mtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, scale, Hm, A, 0.0);
    A = Mtmp;
  }
  Ts += A;

  // Squaring
  for (int s = 0; s < n_sq; ++s) {
    // TODO: Can this work without tmp output storage?
    Mtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ts, Ts, 0.0);
    Ts = Mtmp;
  }

  // unshift and extract first column
  // even for a complex Scalar type, this method produces a real result
  Teuchos::ArrayRCP<Scalar> dd_phi = Teuchos::arcp<Scalar>(lejaOrder);
  const Scalar exp_mu = std::exp(mu); // TODO: is there a Teuchos exp
  for (int i = 0; i < lejaOrder; ++i)
  {
    dd_phi[i] = exp_mu * Ts(i, 0);
  }
  return dd_phi;
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

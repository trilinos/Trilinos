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
      "DD Method to use. 0 for Recurrence. 1 for Taylor Series. 2 for DD_phi");

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
                     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
                     const Scalar cdt
                     )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      phi_order != 0,
      std::invalid_argument,
      "LinOpPhi: phi_order must be zero.");

  const int expansionOrder = this->getExpansionOrder();

  // TODO: optional fractional step size exp(tau*A_tilde)*v
  //const Scalar tau = 1.0;

  // phi_k(L) * v is in range(L)
  const auto rangeSpace = L->range();

  // compute shift and scale parameters
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();
  // TODO: properly apply cdt factor to leja points and shift scale
  // shift *= cdt;
  // scale *= cdt;

  // TODO: update the divided differences (or read from cache)
  //       this should depend on cdt, but that info is not passed down here
  // Get divided differences, we need one more than the polynomial order
  auto lp_dd = getDividedDiffs(phi_order, 1.0, expansionOrder+1);
  TEUCHOS_ASSERT(lp_dd.size() == expansionOrder+1);

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
  auto coeff = lp_dd[0];
  Scalar coeff_re = Scalar(coeff.real());

  //std::cout << "c[0]: " << coeff_re << std::endl;

  Thyra::V_V(vm_k.ptr(), *v);
  Thyra::V_StV(v, coeff_re, *vm_k);

  // storage for error est
  Scalar norm_vm_k = Thyra::norm_inf(*vm_k);
  // norm of the update
  Scalar norm_d_k = std::abs(coeff_re) * norm_vm_k;
  // upper bound on solution size
  Scalar overflow = norm_d_k;
  Thyra::SolveStatus<Scalar> sStatus;

  // leja point index k starts at 1
  int k = 1;
  // leja polynomial term index lp_k starts at 0
  for (int lp_k = 0; k < lp_dd.size() && lp_k < lp_.size(); k++, lp_k++)
  {
    // print the update vector vm_k
    //v->describe(*this->getOStream(), Teuchos::VERB_EXTREME);
    //std::cout << "Norm d_k: " << norm_d_k << " v_k: " << norm_vm_k << std::endl;

    // compute shifted and scaled leja point
    LejaPoint lp_sc = getLpSc(lp_k);
    const Scalar lp_sc_re = Scalar( lp_sc.lp.real() );
    const Scalar lp_sc_im = Scalar( lp_sc.lp.imag() );

    // Real leja point case
    if (lp_sc.lpt == LpType::LPREAL) {
      // extract divided diff
      coeff = lp_dd[k];
      coeff_re = Scalar(coeff.real());
      // std::cout << "c,lp,shift,scale: " << coeff << " " << lp_sc.lp << " real " << shift << " " << scale << std::endl;

      // copy vm_k to temp vector
      Thyra::V_V(av.ptr(), *vm_k);
      // vm_k = (L@vm_k - lp_re[k-1]*vm_k) / scale
      Thyra::apply(*L, Thyra::NOTRANS, *av, vm_k.ptr(), 1/scale, -lp_sc_re/scale);

      // add vm_k*coeff to the final result
      Thyra::Vp_StV(v, coeff_re, *vm_k);

      norm_vm_k = Thyra::norm_inf(*vm_k);
      norm_d_k = std::abs(coeff_re) * norm_vm_k;
      overflow += norm_d_k;
    }
    else if (lp_sc.lpt == LpType::LPCONJ) {
      // extract divided diff
      coeff = lp_dd[k];
      coeff_re = Scalar(coeff.real());
      //std::cout << "c,lp,shift,scale: " << coeff << " " << lp_sc.lp << " " << shift << " " << scale << std::endl;

      // copy vm_k to qm_k vector to save it
      Thyra::V_V(qm_k.ptr(), *vm_k);

      // vm_k = (L@vm_k - lp_re*vm_k) / scale
      Thyra::apply(*L, Thyra::NOTRANS, *qm_k, vm_k.ptr(), 1/scale, -lp_sc_re/scale);

      // add vm_k*coeff to the final result
      Thyra::Vp_StV(v, coeff_re, *vm_k);
      norm_vm_k = Thyra::norm_inf(*vm_k);
      norm_d_k = std::abs(coeff_re) * norm_vm_k;
      overflow += norm_d_k;

      // increment polynomial degree, but keep Leja point and handle conjugate pair
      k++;
      // std::cout << "qm" << std::endl;
      // vm_k->describe(*this->getOStream(), Teuchos::VERB_EXTREME);

      if (k < lp_dd.size())
      {
        // conjugate update
        coeff = lp_dd[k];
        coeff_re = Scalar(coeff.real());

        //std::cout << "Norm d_k: " << norm_d_k << " v_k: " << norm_vm_k << std::endl;
        //std::cout << "c,lp: " << coeff << " " << std::conj(lp_sc.lp) << std::endl;

        // copy vm_k to a new temp vector (don't overwrite qm_k)
        Thyra::V_V(av.ptr(), *vm_k);
        // vm_k = (L@vm_k - lp_re*vm_k) / scale + ((lp_im/scale)**2)*qm_k
        Thyra::apply(*L, Thyra::NOTRANS, *av, vm_k.ptr(), 1/scale, -lp_sc_re/scale);
        Thyra::Vp_StV(vm_k.ptr(), (lp_sc_im/scale) * (lp_sc_im/scale), *qm_k);

        // add vm_k*coeff to the final result
        Thyra::Vp_StV(v, coeff_re, *vm_k);
        norm_vm_k = Thyra::norm_inf(*vm_k);
        norm_d_k = std::abs(coeff_re) * norm_vm_k;
        overflow += norm_d_k;
      }
    }
    else {
      // TODO: ERROR
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

    // terminate if the update drops below likely significance
    //if (norm_vm_k < overflow / this->leja_tol_)
    // terminate if the update drops below user tol
    if (k >= expansionOrder || norm_d_k < this->leja_tol_)
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

  std::cout << sStatus.message << std::endl;

  return sStatus;
}


template <class Scalar>
void PhiEvaluatorLeja<Scalar>::initLejaPointsBase()
{
  const int exp_order = this->getExpansionOrder();

  // The maximum number of Leja points is: exp_order + 1, but due to conjugacy, we technically need less points stored:
  // For real Leja points given as the real part of the points below,
  //    this is perfect, since every conjugate Leja point maps to a single real point
  // For imaginary Leja points given as the imaginary part,
  //    this mostly correct but sufficient for exp_order + 1 > 2, since the first two points map both to 0
  // For conjugate Leja points, we only need 2 + exp_order / 2,
  //    this still works since it is an upper bound
  const int maxLejaOrder = std::max(2, exp_order + 1);

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
std::tuple<Scalar, Scalar> PhiEvaluatorLeja<Scalar>::getShiftScale()
{
  // real half axis
  Scalar hx_re = (leja_b_ - leja_a_) / 2.0;
  // imaj half axis
  Scalar hx_im = leja_c_;
  // leja ellipse shift and scale parameters
  Scalar shift = (leja_a_ + leja_b_) / 2.0;
  Scalar scale = (hx_re + hx_im) / 2.0;
  return std::make_tuple(shift, scale);
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setExpansionOrder(int order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      order <= 0,
      std::invalid_argument,
      "setExpansionOrder: order must be positive.");

  expansionOrder_ = order;

  initLejaPointsBase();
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setLejaEllipse(Scalar a, Scalar b, Scalar c)
{
  TEUCHOS_ASSERT(a <= b);
  TEUCHOS_ASSERT(c >= 0.0);
  leja_a_ = a;
  leja_b_ = b;
  leja_c_ = c;

  // update the leja points
  lp_ = Teuchos::arcp<LejaPoint>(lejaPointsBase_.size());
  Scalar hx_re = (leja_b_ - leja_a_) / 2.0;
  Scalar hx_im = leja_c_;
  Scalar scale = (hx_re + hx_im) / 2.0;
  for (int i=0; i < lejaPointsBase_.size(); ++i) {
    auto lp_real = lejaPointsBase_[i].lp.real();
    auto lp_imag = lejaPointsBase_[i].lp.imag();
    auto lp = std::complex(lp_real * hx_re / scale, lp_imag * hx_im / scale);
    lp_[i] = {lp, lejaPointsBase_[i].lpt};
  }
}

template <class Scalar>
LejaPoint PhiEvaluatorLeja<Scalar>::getLpSc(int i)
{
  TEUCHOS_ASSERT(i < lp_.size());
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();
  LejaPoint lp = this->lp_[i];
  LejaPoint lp_sc = LejaPoint{shift + scale * lp.lp, lp.lpt};
  return lp_sc;
}

template <class Scalar>
void PhiEvaluatorLeja<Scalar>::setPhiEvaluatorValues(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  PhiEvaluator<Scalar>::setPhiEvaluatorValues(pl);

  //pl->validateParametersAndSetDefaults(*getValidParameters());

  leja_tol_ = pl->get<double>("leja_tol", 1.0e-18);
  ddMethod_ = pl->get<int>("Leja DD Method", 1);
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
Teuchos::ArrayRCP<std::complex<double>> PhiEvaluatorLeja<Scalar>::getDividedDiffs(const int k, const Scalar cdt, const int exp_order)
{
  if (ddMethod_ == 2) {
    std::cout << "Calling dd_phi divided difference method" << std::endl;
    return getDividedDiffsPhi(k, cdt, exp_order);
  }
  // Default: Taylor series (real arithmetic)
  std::cout << "Calling taylor series divided difference method" << std::endl;
  return getDividedDiffsTSR(k, cdt, exp_order);
}

template <class Scalar>
Teuchos::ArrayRCP<std::complex<double>> PhiEvaluatorLeja<Scalar>::getDividedDiffsPhi(
    const int phi_order, const Scalar cdt, const int exp_order)
{
  // Compute Newton divided differences of phi_{phi_order} at the Leja points using
  // the Zivcovich (2019) H-factorization + scaling-and-squaring method.
  //
  // Ref: F. Zivcovich. "Fast and accurate computation of divided differences for
  //      analytic functions, with an application to the exponential function."
  //      Dolomites Research Notes on Approximation. 12. 2019.
  using cplx = std::complex<double>;
  using cplxl = std::complex<long double>;

  const int n_leja = exp_order;                   // number of output coefficients
  const int l      = phi_order;                   // leading zeros for phi_l
  const int total  = l + n_leja;                  // total interpolation points
  const int n      = total - 1;                   // highest index (0-based)
  const int p      = 30;                          // Taylor truncation terms
  const int cap_n  = n + p;                       // Taylor degree N in the paper

  // get spectrum shift/scale
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();

  // build z[0..total]
  // z[0..l]         = 0  (leading zeros encode the phi_k recurrence)
  // z[l..total]     = cdt * lp_sc(i), expanding LPCONJ pairs to (lp, conj(lp))
  std::vector<cplxl> z(total, cplx(0.0, 0.0));
  {
    int z_idx = l;
    for (int lp_idx = 0; lp_idx < lp_.size() && z_idx < total; ++lp_idx)
    {
      LejaPoint lp_sc = getLpSc(lp_idx);
      if (lp_sc.lpt == LPCONJ)
      {
        // upper-half-plane point
        z[z_idx++] = cplx(cdt * lp_sc.lp.real(),  cdt * lp_sc.lp.imag());
        // conjugate (lower-half-plane)
        if (z_idx < total)
          z[z_idx++] = cplx(cdt * lp_sc.lp.real(), -cdt * lp_sc.lp.imag());
      }
      else  // LPREAL: zero out any floating-point imaginary noise
      {
        z[z_idx++] = cplx(cdt * lp_sc.lp.real(), 0.0);
      }
    }
  }

  // shift by mean mu
  cplx mu(0.0, 0.0);
  for (int i = 0; i < total; ++i) mu += z[i];
  mu /= double(total);
  for (int i = 0; i < total; ++i) z[i] -= mu;

  // build lower-triangle of F
  Teuchos::SerialDenseMatrix<int, cplx> F_mat(total, total);
  F_mat.putScalar(cplx(0.0, 0.0));
  for (int i0 = 0; i0 < n; ++i0)
    for (int j0 = i0 + 1; j0 <= n; ++j0)
      F_mat(j0, i0) = z[i0] - z[j0];

  // number of squarings
  double max_abs = 0.0;
  for (int i0 = 0; i0 < n; ++i0)
    for (int j0 = i0 + 1; j0 <= n; ++j0)
    {
      double v = std::abs(F_mat(j0, i0));
      if (v > max_abs) max_abs = v;
    }
  const int    s     = std::max(int(std::ceil(max_abs / 3.5)), 1);
  const double s_dbl = double(s);

  // seed dd[0..cap_n]: dd[kk] = 1 / (kk! * s^kk)
  std::vector<cplxl> dd(cap_n + 1, cplxl(0.0, 0.0));
  dd[0] = cplx(1.0, 0.0);
  // avoid overflow
  long double running_denom = 1.0;
  for (int kk = 1; kk <= cap_n; ++kk)
  {
    running_denom *= kk * s_dbl;
    dd[kk] = cplxl(1.0 / running_denom, 0.0);
  }

  // H-factorization sweep
  for (int j = n; j >= 0; --j)
  {
    // First inner loop: Taylor remainder sweep
    // k0 from cap_n-1 downto n-j (inclusive): dd[k0] += z[j] * dd[k0+1]
    for (int k0 = cap_n - 1; k0 >= (n - j); --k0)
      dd[k0] = dd[k0] + z[j] * dd[k0 + 1];

    // Second inner loop: divided-difference sweep using F lower triangle
    // k0 from (n-j-1) downto 0: dd[k0] += F(k0+j+1, j) * dd[k0+1]
    for (int k0 = (n - j) - 1; k0 >= 0; --k0)
      dd[k0] = dd[k0] + (cplxl) F_mat(k0 + j + 1, j) * dd[k0 + 1];

    // Store dd[0..=n-j] into upper-triangle row j of F
    for (int col = 0; col <= (n - j); ++col)
      F_mat(j, j + col) = dd[col];
  }

  // overwrite diagonal: F(i,i) = exp(z[i] / s)
  for (int i = 0; i <= n; ++i)
    F_mat(i, i) = std::exp(cplx(z[i].real() / s_dbl, z[i].imag() / s_dbl));

  // zero lower triangle
  for (int i = 1; i <= n; ++i)
    for (int j = 0; j < i; ++j)
      F_mat(i, j) = cplx(0.0, 0.0);

  // squaring
  Teuchos::SerialDenseMatrix<int, cplx> dd_row(1, total);
  for (int j = 0; j < total; ++j)
    dd_row(0, j) = F_mat(0, j);

  Teuchos::SerialDenseMatrix<int, cplx> tmp_row(1, total);
  for (int sq = 0; sq < s - 1; ++sq)
  {
    tmp_row.putScalar(cplx(0.0, 0.0));
    // tmp_row = 1.0 * dd_row * F_mat + 0.0 * tmp_row
    tmp_row.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, cplx(1.0, 0.0), dd_row, F_mat, cplx(0.0, 0.0));
    dd_row = tmp_row;
  }

  // out[i] = exp(mu) * scale^i * dd_row[l + i]
  Teuchos::ArrayRCP<cplx> out = Teuchos::arcp<cplx>(n_leja);
  const cplx exp_mu   = std::exp(mu);
  const cplx scale_c  = cplx(double(scale), 0.0);
  cplx       scale_pw = cplx(1.0, 0.0);
  for (int i = 0; i < n_leja; ++i)
  {
    out[i]   = exp_mu * scale_pw * dd_row(0, l + i);
    scale_pw *= scale_c;
  }

  return out;
}


template <class Scalar>
Teuchos::ArrayRCP<std::complex<double>> PhiEvaluatorLeja<Scalar>::getDividedDiffsRC(const int phi_order, const Scalar cdt, const int exp_order)
{
  TEUCHOS_ASSERT(phi_order == 0);

  //TODO: as long as leja points use double, shift and scale should too
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();
  auto transform = [&](std::complex<double> x) {
    return cdt * (shift + scale * x);
  };

  Teuchos::ArrayRCP<std::complex<double>> x = Teuchos::arcp<std::complex<double>>(exp_order);
  Teuchos::ArrayRCP<std::complex<double>> d_x = Teuchos::arcp<std::complex<double>>(exp_order);

  // initialize list of Leja points and function values
  for (int idx = 0, lp_idx = 0; idx < exp_order && lp_idx < lp_.size(); idx++, lp_idx++)
  {
    LejaPoint lp = this->lp_[lp_idx];
    if (lp.lpt == LPCONJ)
    {
      x[idx] = lp.lp;
      d_x[idx] = std::exp(transform(x[idx]));
      if (++idx < exp_order)
      {
        x[idx] = std::conj(lp.lp);
        d_x[idx] = std::exp(transform(x[idx]));
      }
    }
    else
    {
      x[idx] = lp.lp;
      d_x[idx] = std::exp(transform(x[idx]));
    }
  }

  for (int idx = 0; idx < exp_order-1; idx++)
  {
    // Compute the next set of divided differences
    for (int idy = idx+1; idy < exp_order; idy++)
    {
      d_x[idy] = (d_x[idy] - d_x[idx]) / (x[idy] - x[idx]);
    }
  }

  return d_x;
}

template <class Scalar>
Teuchos::ArrayRCP<std::complex<double>> PhiEvaluatorLeja<Scalar>::getDividedDiffsTS(const int phi_order, const Scalar cdt, const int exp_order)
{
  TEUCHOS_ASSERT(phi_order == 0);

  int m = exp_order;
  Teuchos::ArrayRCP<std::complex<double>> out = Teuchos::arcp<std::complex<double>>(m);

#ifdef HAVE_TEUCHOS_COMPLEX
  // get shift and scale parameters
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();

  // build the shifted and scaled Hm matrix
  Teuchos::SerialDenseMatrix<int, std::complex<double>> Hm(m, m);
  // diagonal elements are the leja points
  int dd_idx = 0;
  int lp_idx = 0;
  while (dd_idx < m) {
    LejaPoint lp_sc = getLpSc(lp_idx);
    // conj lp case
    if (lp_sc.lpt == LPCONJ) {
      if (dd_idx == m) break;
      Hm(dd_idx, dd_idx) = cdt*lp_sc.lp;
      if (dd_idx+1 < m) Hm(dd_idx+1, dd_idx) = cdt*scale;
      dd_idx += 1;
      if (dd_idx == m) break;
      Hm(dd_idx, dd_idx) = cdt*std::conj(lp_sc.lp);
      if (dd_idx+1 < m) Hm(dd_idx+1, dd_idx) = cdt*scale;
      dd_idx += 1;
    }
    else {
      if (dd_idx == m) break;
      Hm(dd_idx, dd_idx) = cdt*lp_sc.lp;
      if (dd_idx+1 < m) Hm(dd_idx+1, dd_idx) = cdt*scale;
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
  std::complex<double> diag_sum = std::complex(0.0, 0.0);
  for (int i=0; i < m; ++i) {
    diag_sum += Hm(i, i);
  }
  std::complex<double> mu = diag_sum / double(m);

  // shift diagonal to zero mean
  for (int i=0; i < m; ++i) {
    Hm(i, i) -= mu;
  }

  // Scaling
  double s_scale = Hm.normInf();
  int n_sq = std::max(int( std::ceil((std::log(s_scale) - std::log(2.0)) / std::log(2.0)) ), 1);

  //n_sq += 2; // increase number of scalings to reduce Taylor poly size.

  double h_scale = 1.0 / std::pow(2.0, n_sq);
  Hm.scale(h_scale);

  // compute phi_0(Hm) by Taylor series
  //copy Hm to A
  Teuchos::SerialDenseMatrix<int, std::complex<double>> A(Teuchos::Copy, Hm);

  Teuchos::SerialDenseMatrix<int, std::complex<double>> Ts(m, m);
  Ts = 0.;

  //auto fact = [](double n) -> double { return std::tgamma(1.0 + n); };

  for (int i=0; i < m; ++i) {
    Ts(i, i) = std::complex(1.0, 0.0);
  }

  // Ts = I/(p!) + Hm^1/(1+p)! + Hm^2/(2+p)! ...
  int ts_order = 17;

  //std::cout << "ts_order: " << ts_order << std::endl;
  //std::cout << "s_scale: " << s_scale << std::endl;
  //std::cout << "mu: " << mu << std::endl;
  //std::cout << "n_sq: " << n_sq << std::endl;
  //std::cout << "h_scale: " << h_scale << std::endl;

  Teuchos::SerialDenseMatrix<int, std::complex<double>> Mtmp(m, m);
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
  // std::cout << "Ts sq(0, 0): " << Ts(0, 0) << std::endl;
  // std::cout << "Ts sq(1, 0): " << Ts(1, 0) << std::endl;
  // std::cout << "Ts sq(2, 0): " << Ts(2, 0) << std::endl;
  // std::cout << "Ts sq(3, 0): " << Ts(3, 0) << std::endl;

  // unshift and extract first column
  for (int i=0; i < m; ++i) {
    out[i] = std::exp(mu) * Ts(i, 0);
  }
#else
  std::cout << "WARNING: getDividedDiffsTS requires Trilinos configured with Trilinos_ENABLE_COMPLEX=ON" << std::endl;
  std::cout << "WARNING: falling back to getDividedDiffs." << std::endl;
  // TODO: implement fallback dd implementation here.
  return getDividedDiffsRC(phi_order, cdt);
#endif
  return out;
}

template <class Scalar>
Teuchos::ArrayRCP<std::complex<double>> PhiEvaluatorLeja<Scalar>::getDividedDiffsTSR(const int phi_order, const Scalar cdt, const int exp_order)
{
  TEUCHOS_ASSERT(phi_order == 0);

  const int m = exp_order;

  //TODO: still need a complex output for compatibility although result is real
  Teuchos::ArrayRCP<std::complex<double>> out = Teuchos::arcp<std::complex<double>>(m);

  // get shift and scale parameters
  Scalar shift, scale;
  std::tie(shift, scale) = getShiftScale();

  // build the shifted and scaled Hm matrix
  Teuchos::SerialDenseMatrix<int, Scalar> Hm(m, m);
  // diagonal elements are the leja points

  for (int lp_idx = 0, dd_idx = 0; lp_idx < lp_.size() && dd_idx < m; lp_idx++, dd_idx++) {
    LejaPoint lp_sc = getLpSc(lp_idx);
    // conj lp case
    if (lp_sc.lpt == LPCONJ) {
      Hm(dd_idx, dd_idx) = cdt * lp_sc.lp.real();
      if (dd_idx + 1 < m) Hm(dd_idx + 1, dd_idx) = cdt * scale;

      if (++dd_idx < m) {
        Hm(dd_idx - 1, dd_idx) = - cdt * lp_sc.lp.imag() * lp_sc.lp.imag() / scale;
        Hm(dd_idx, dd_idx) = cdt * lp_sc.lp.real();
        if (dd_idx + 1 < m)
          Hm(dd_idx + 1, dd_idx) = cdt * scale;
      }
    }
    else {
      Hm(dd_idx, dd_idx) = cdt * lp_sc.lp.real();
      if (dd_idx + 1 < m)
        Hm(dd_idx + 1, dd_idx) = cdt * scale;
    }
  }

  // compute diagonal mean
  Scalar diag_sum = 0;
  for (int i = 0; i < m; ++i) {
    diag_sum += Hm(i, i);
  }
  double mu = diag_sum / double(m);

  // shift diagonal to zero mean
  for (int i = 0; i < m; ++i) {
    Hm(i, i) -= mu;
  }

  // Scaling
  double s_scale = Hm.normInf();
  int n_sq = std::max(int( std::ceil((std::log(s_scale) - std::log(2.0)) / std::log(2.0)) ), 1);

  //n_sq += 2; // increase number of scalings to reduce Taylor poly size.

  double h_scale = 1.0 / std::pow(2.0, n_sq);
  Hm.scale(h_scale);

  // compute phi_0(Hm) by Taylor series
  //copy Hm to A
  Teuchos::SerialDenseMatrix<int, Scalar> A(Teuchos::Copy, Hm);

  Teuchos::SerialDenseMatrix<int, Scalar> Ts(m, m);
  Ts = 0.;

  //auto fact = [](double n) -> double { return std::tgamma(1.0 + n); };

  for (int i = 0; i < m; ++i) {
    Ts(i, i) = 1.0;
  }

  // Ts = I/(p!) + Hm^1/(1+p)! + Hm^2/(2+p)! ...
  int ts_order = 17;

  //std::cout << "ts_order: " << ts_order << std::endl;
  //std::cout << "s_scale: " << s_scale << std::endl;
  //std::cout << "mu: " << mu << std::endl;
  //std::cout << "n_sq: " << n_sq << std::endl;
  //std::cout << "h_scale: " << h_scale << std::endl;

  Teuchos::SerialDenseMatrix<int, Scalar> Mtmp(m, m);
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
  for (int s = 0; s < n_sq; ++s) {
    // TODO: Can this work without tmp output storage?
    Mtmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ts, Ts, 0.0);
    Ts = Mtmp;
  }

  // unshift and extract first column
  for (int i = 0; i < m; ++i) {
    out[i] = std::exp(mu) * Ts(i, 0);
  }

  return out;
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

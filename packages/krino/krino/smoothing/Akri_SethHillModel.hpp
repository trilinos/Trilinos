// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_SETHHILLMODEL_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_SETHHILLMODEL_HPP_

#include <limits>
#include <utility>

#include <Akri_KinematicUtils.hpp>
#include <Akri_SethHillConfig.hpp>
#include <MiniTensor.h>

namespace krino {

// Ported from Norma.jl (https://github.com/sandialabs/Norma.jl)
// by Alejandro Mota

namespace sethHillUtils {

// Matrix power A^p for integer p.
// Uses binary_powering for |p| and inverse for negative p.
template<minitensor::Index N>
minitensor::Tensor<double, N> power(const minitensor::Tensor<double, N> & A, int p)
{
  if (p == 0)
    return minitensor::identity<double, N>();
  if (p == 1)
    return A;
  if (p == -1)
    return minitensor::inverse(A);
  if (p > 0)
    return minitensor::binary_powering(A, static_cast<minitensor::Index>(p));
  // p < -1: inverse of A^|p|
  return minitensor::inverse(
      minitensor::binary_powering(A, static_cast<minitensor::Index>(-p)));
}

// Matrix power A^p for integer p.
// Uses binary_powering for |p| and inverse for negative p.
template<minitensor::Index N>
std::pair<bool,minitensor::Tensor<double, N>> safe_power(const minitensor::Tensor<double, N> & A, int p)
{
  if (p == 0)
    return std::make_pair(true, minitensor::identity<double, N>());
  if (p == 1)
    return std::make_pair(true,A);
  if (p == -1)
  {
    return kinematicUtils::safe_inverse(A);
  }
  if (p > 0)
    return std::make_pair(true, minitensor::binary_powering(A, static_cast<minitensor::Index>(p)));
  // p < -1: inverse of A^|p|
  return kinematicUtils::safe_inverse(minitensor::binary_powering(A, static_cast<minitensor::Index>(-p)));
}

}

// Seth-Hill energy density W(F, J) and its first Piola-Kirchhoff stress P = dW/dF, in 3D and 2D.
class SethHillModel
{
public:
  using Tensor = minitensor::Tensor<double, 3>;
  using Tensor2d = minitensor::Tensor<double, 2>;

  explicit SethHillModel(const SethHillParams & params) : myParams(params) {}

  double energy(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const double Jm23 = kinematicUtils::pow_m23(J);
    const Tensor Cbar = Jm23 * C;

    const double Jm  = std::pow(J, myParams.m);
    const double Jnm = 1.0 / Jm;

    const auto [successCbar_n, Cbar_n] = sethHillUtils::safe_power(Cbar, myParams.n);
    const auto [successCbar_nn, Cbar_nn] = sethHillUtils::safe_power(Cbar, -myParams.n);
    if (!successCbar_n || !successCbar_nn)
      return std::numeric_limits<double>::max();

    const Tensor Cbar_2n  = Cbar_n * Cbar_n;
    const Tensor Cbar_n2n = Cbar_nn * Cbar_nn;

    const double trCbar_n   = minitensor::trace(Cbar_n);
    const double trCbar_nn  = minitensor::trace(Cbar_nn);
    const double trCbar_2n  = minitensor::trace(Cbar_2n);
    const double trCbar_n2n = minitensor::trace(Cbar_n2n);

    const double m2 = static_cast<double>(myParams.m * myParams.m);
    const double n2 = static_cast<double>(myParams.n * myParams.n);

    const double Wbulk  = myParams.kappa / (4.0 * m2) *
        ((Jm - 1.0) * (Jm - 1.0) + (Jnm - 1.0) * (Jnm - 1.0));
    const double Wshear = myParams.mu / (4.0 * n2) *
        (trCbar_2n + trCbar_n2n - 2.0 * trCbar_n - 2.0 * trCbar_nn + 6.0);

    return Wbulk + Wshear;
  }

  Tensor stress(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const Tensor Finv = minitensor::inverse(F);
    const Tensor FinvT = minitensor::transpose(Finv);

    const double Jm   = std::pow(J, myParams.m);
    const double Jnm  = 1.0 / Jm;
    const double J2m  = Jm * Jm;
    const double Jn2m = 1.0 / J2m;

    const double Jm23 = kinematicUtils::pow_m23(J);
    const double Jp23 = kinematicUtils::pow_p23(J);
    const Tensor Cbar     = Jm23 * C;
    const Tensor Cbar_inv = Jp23 * (Finv * FinvT);

    const Tensor Cbar_n   = sethHillUtils::power(Cbar, myParams.n);
    const Tensor Cbar_nn  = sethHillUtils::power(Cbar_inv, myParams.n);
    const Tensor Cbar_2n  = Cbar_n * Cbar_n;
    const Tensor Cbar_n2n = Cbar_nn * Cbar_nn;

    const double trCbar_n   = minitensor::trace(Cbar_n);
    const double trCbar_nn  = minitensor::trace(Cbar_nn);
    const double trCbar_2n  = minitensor::trace(Cbar_2n);
    const double trCbar_n2n = minitensor::trace(Cbar_n2n);

    // Pbulk = kappa/(2m) * (J^2m - J^m - J^-2m + J^-m) * F^-T
    const double bulk_coeff = myParams.kappa / (2.0 * myParams.m) *
        (J2m - Jm - Jn2m + Jnm);
    const Tensor Pbulk = bulk_coeff * FinvT;

    // Pshear = mu/n * [1/3*(-trCbar2n + trCbarn + trCbar-2n - trCbar-n)*F^-T
    //                   + F^-T * (Cbar2n - Cbarn - Cbar-2n + Cbar-n)]
    const double iso_coeff = (1.0 / 3.0) *
        (-trCbar_2n + trCbar_n + trCbar_n2n - trCbar_nn);
    const Tensor dev_part = (Cbar_2n - Cbar_n) - (Cbar_n2n - Cbar_nn);
    const Tensor Pshear = (myParams.mu / myParams.n) *
        (iso_coeff * FinvT + FinvT * dev_part);

    return Pbulk + Pshear;
  }

  double energy(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const double Jinv = 1.0 / J;
    const Tensor2d Cbar = Jinv * C;

    const double Jm  = std::pow(J, myParams.m);
    const double Jnm = 1.0 / Jm;

    const auto [successCbar_n,  Cbar_n]  = sethHillUtils::safe_power(Cbar,  myParams.n);
    const auto [successCbar_nn, Cbar_nn] = sethHillUtils::safe_power(Cbar, -myParams.n);
    if (!successCbar_n || !successCbar_nn)
      return std::numeric_limits<double>::max();

    const Tensor2d Cbar_2n  = Cbar_n * Cbar_n;
    const Tensor2d Cbar_n2n = Cbar_nn * Cbar_nn;

    const double trCbar_n   = minitensor::trace(Cbar_n);
    const double trCbar_nn  = minitensor::trace(Cbar_nn);
    const double trCbar_2n  = minitensor::trace(Cbar_2n);
    const double trCbar_n2n = minitensor::trace(Cbar_n2n);

    const double m2 = static_cast<double>(myParams.m * myParams.m);
    const double n2 = static_cast<double>(myParams.n * myParams.n);

    const double Wbulk  = myParams.kappa / (4.0 * m2) *
        ((Jm - 1.0) * (Jm - 1.0) + (Jnm - 1.0) * (Jnm - 1.0));

    const double Wshear = myParams.mu / (4.0 * n2) *
        (trCbar_2n + trCbar_n2n - 2.0 * trCbar_n - 2.0 * trCbar_nn + 4.0);

    return Wbulk + Wshear;
  }

  Tensor2d stress(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const Tensor2d Finv = minitensor::inverse(F);
    const Tensor2d FinvT = minitensor::transpose(Finv);

    const double Jm   = std::pow(J, myParams.m);
    const double Jnm  = 1.0 / Jm;
    const double J2m  = Jm * Jm;
    const double Jn2m = 1.0 / J2m;

    const double Jinv = 1.0 / J;
    const Tensor2d Cbar     = Jinv * C;
    const Tensor2d Cbar_inv = J * (Finv * FinvT);

    const Tensor2d Cbar_n   = sethHillUtils::power(Cbar, myParams.n);
    const Tensor2d Cbar_nn  = sethHillUtils::power(Cbar_inv, myParams.n);
    const Tensor2d Cbar_2n  = Cbar_n * Cbar_n;
    const Tensor2d Cbar_n2n = Cbar_nn * Cbar_nn;

    const double trCbar_n   = minitensor::trace(Cbar_n);
    const double trCbar_nn  = minitensor::trace(Cbar_nn);
    const double trCbar_2n  = minitensor::trace(Cbar_2n);
    const double trCbar_n2n = minitensor::trace(Cbar_n2n);

    // Pbulk = kappa/(2m) * (J^2m - J^m - J^-2m + J^-m) * F^-T
    const double bulk_coeff = myParams.kappa / (2.0 * myParams.m) *
        (J2m - Jm - Jn2m + Jnm);
    const Tensor2d Pbulk = bulk_coeff * FinvT;

    // 2D:
    // Pshear = mu/n * [1/2*(-trCbar2n + trCbarn + trCbar-2n - trCbar-n)*F^-T
    //                  + F^-T * (Cbar2n - Cbarn - Cbar-2n + Cbar-n)]
    const double iso_coeff = 0.5 *
        (-trCbar_2n + trCbar_n + trCbar_n2n - trCbar_nn);
    const Tensor2d dev_part = (Cbar_2n - Cbar_n) - (Cbar_n2n - Cbar_nn);
    const Tensor2d Pshear = (myParams.mu / myParams.n) *
        (iso_coeff * FinvT + FinvT * dev_part);

    return Pbulk + Pshear;
  }

private:
  SethHillParams myParams;
};

// Seth-Hill specialized and optimized for m = 1, n = 1, kappa = 1, mu = 1.
class SimplifiedSethHillModel
{
public:
  using Tensor = minitensor::Tensor<double, 3>;
  using Tensor2d = minitensor::Tensor<double, 2>;

  // Strain energy W(F) specialized and optimized for
  double energy(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;

    const Tensor Finv  = minitensor::inverse(F);
    const Tensor FinvT = minitensor::transpose(Finv);

    // C^{-1} = F^{-1} F^{-T}
    const Tensor Cinv = Finv * FinvT;

    const Tensor C2    = C * C;
    const Tensor Cinv2 = Cinv * Cinv;

    const double trC     = minitensor::trace(C);
    const double trC2    = minitensor::trace(C2);
    const double trCinv  = minitensor::trace(Cinv);
    const double trCinv2 = minitensor::trace(Cinv2);

    const double Jinv = 1.0 / J;

    const double Jm23 = kinematicUtils::pow_m23(J);
    const double Jp23 = kinematicUtils::pow_p23(J);

    const double Jm43 = Jm23 * Jm23;
    const double Jp43 = Jp23 * Jp23;

    // Wbulk = 1/4 * (J - 1)^2 * (1 + J^{-2})
    const double Jminus1 = J - 1.0;
    const double Wbulk = 0.25 * Jminus1 * Jminus1 * (1.0 + Jinv * Jinv);

    // Wshear = 1/4 * [
    //   J^{-4/3} tr(C^2)
    // + J^{ 4/3} tr(C^{-2})
    // - 2 J^{-2/3} tr(C)
    // - 2 J^{ 2/3} tr(C^{-1})
    // + 6
    // ]
    const double Wshear = 0.25 *
        (Jm43 * trC2 + Jp43 * trCinv2
         - 2.0 * Jm23 * trC
         - 2.0 * Jp23 * trCinv
         + 6.0);

    return Wbulk + Wshear;
  }

  // 1st Piola-Kirchhoff stress P = dW/dF
  Tensor stress(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;

    const Tensor Finv  = minitensor::inverse(F);
    const Tensor FinvT = minitensor::transpose(Finv);

    // C^{-1} = F^{-1} F^{-T}
    const Tensor Cinv = Finv * FinvT;

    const Tensor C2    = C * C;
    const Tensor Cinv2 = Cinv * Cinv;

    const double trC     = minitensor::trace(C);
    const double trC2    = minitensor::trace(C2);
    const double trCinv  = minitensor::trace(Cinv);
    const double trCinv2 = minitensor::trace(Cinv2);

    const double Jinv = 1.0 / J;

    const double Jm23 = kinematicUtils::pow_m23(J);
    const double Jp23 = kinematicUtils::pow_p23(J);

    const double Jm43 = Jm23 * Jm23;
    const double Jp43 = Jp23 * Jp23;

    // Pbulk = 1/2 * (J - 1) * (J + J^{-2}) F^{-T}
    const double bulk_coeff =
        0.5 * (J - 1.0) * (J + Jinv * Jinv);

    const Tensor Pbulk = bulk_coeff * FinvT;

    // The deviatoric tensor before left multiplication by F^{-T} is:
    //
    // dev_part =
    //   J^{-4/3} C^2
    // - J^{-2/3} C
    // - J^{ 4/3} C^{-2}
    // + J^{ 2/3} C^{-1}
    //
    // iso_coeff = -1/3 tr(dev_part)
    const double tr_dev =
        Jm43 * trC2
        - Jm23 * trC
        - Jp43 * trCinv2
        + Jp23 * trCinv;

    const double iso_coeff = -(1.0 / 3.0) * tr_dev;

    // Use identities:
    // F^{-T} C  = F
    // F^{-T} C^2 = F C
    //
    // Therefore:
    // F^{-T} dev_part =
    //   J^{-4/3} F C
    // - J^{-2/3} F
    // - J^{ 4/3} F^{-T} C^{-2}
    // + J^{ 2/3} F^{-T} C^{-1}
    const Tensor Pshear =
        iso_coeff * FinvT
        + Jm43 * (F * C)
        - Jm23 * F
        - Jp43 * (FinvT * Cinv2)
        + Jp23 * (FinvT * Cinv);

    return Pbulk + Pshear;
  }

  // 2D strain energy W(F) specialized and optimized for
  // m = 1, n = 1, kappa = 1, mu = 1
  //
  double energy(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;

    const Tensor2d Finv  = minitensor::inverse(F);
    const Tensor2d FinvT = minitensor::transpose(Finv);

    // C^{-1} = F^{-1} F^{-T}
    const Tensor2d Cinv = Finv * FinvT;

    const Tensor2d C2    = C * C;
    const Tensor2d Cinv2 = Cinv * Cinv;

    const double trC     = minitensor::trace(C);
    const double trC2    = minitensor::trace(C2);
    const double trCinv  = minitensor::trace(Cinv);
    const double trCinv2 = minitensor::trace(Cinv2);

    const double Jinv  = 1.0 / J;
    const double Jinv2 = Jinv * Jinv;
    const double J2    = J * J;

    // Wbulk = 1/4 * [ (J - 1)^2 + (J^{-1} - 1)^2 ]
    //       = 1/4 * (J - 1)^2 * (1 + J^{-2})
    const double Jminus1 = J - 1.0;
    const double Wbulk = 0.25 * Jminus1 * Jminus1 * (1.0 + Jinv2);

    // In 2D:
    // Cbar      = J^{-1} C
    // Cbar^{-1} = J C^{-1}
    //
    // Wshear = 1/4 * [
    //   tr(Cbar^2)
    // + tr(Cbar^{-2})
    // - 2 tr(Cbar)
    // - 2 tr(Cbar^{-1})
    // + 4
    // ]
    //
    // = 1/4 * [
    //   J^{-2} tr(C^2)
    // + J^2 tr(C^{-2})
    // - 2 J^{-1} tr(C)
    // - 2 J tr(C^{-1})
    // + 4
    // ]
    const double Wshear = 0.25 *
        (Jinv2 * trC2
         + J2 * trCinv2
         - 2.0 * Jinv * trC
         - 2.0 * J * trCinv
         + 4.0);

    return Wbulk + Wshear;
  }

  // 2D 1st Piola-Kirchhoff stress P = dW/dF
  // specialized and optimized for
  // m = 1, n = 1, kappa = 1, mu = 1
  //
  Tensor2d stress(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;

    const Tensor2d Finv  = minitensor::inverse(F);
    const Tensor2d FinvT = minitensor::transpose(Finv);

    // C^{-1} = F^{-1} F^{-T}
    const Tensor2d Cinv = Finv * FinvT;

    const Tensor2d C2    = C * C;
    const Tensor2d Cinv2 = Cinv * Cinv;

    const double trC     = minitensor::trace(C);
    const double trC2    = minitensor::trace(C2);
    const double trCinv  = minitensor::trace(Cinv);
    const double trCinv2 = minitensor::trace(Cinv2);

    const double Jinv  = 1.0 / J;
    const double Jinv2 = Jinv * Jinv;
    const double J2    = J * J;

    // Pbulk = 1/2 * (J^2 - J - J^{-2} + J^{-1}) F^{-T}
    //
    // Equivalent factored form:
    //
    // Pbulk = 1/2 * (J - 1) * (J + J^{-2}) F^{-T}
    const double bulk_coeff = 0.5 * (J - 1.0) * (J + Jinv2);
    const Tensor2d Pbulk = bulk_coeff * FinvT;

    // In 2D:
    //
    // dev_part = Cbar^2 - Cbar - Cbar^{-2} + Cbar^{-1}
    // = J^{-2} C^2
    // - J^{-1} C
    // - J^2 C^{-2}
    // + J C^{-1}
    //
    // iso_coeff = -1/2 * tr(dev_part)
    const double tr_dev =
        Jinv2 * trC2
        - Jinv * trC
        - J2 * trCinv2
        + J * trCinv;

    const double iso_coeff = -0.5 * tr_dev;

    // Use identities:
    //
    // F^{-T} C   = F
    // F^{-T} C^2 = F C
    //
    // Therefore:
    //
    // F^{-T} dev_part =
    //   J^{-2} F C
    // - J^{-1} F
    // - J^2 F^{-T} C^{-2}
    // + J F^{-T} C^{-1}
    const Tensor2d Pshear =
        iso_coeff * FinvT
        + Jinv2 * (F * C)
        - Jinv * F
        - J2 * (FinvT * Cinv2)
        + J * (FinvT * Cinv);

    return Pbulk + Pshear;
  }
};

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_SETHHILLMODEL_HPP_ */

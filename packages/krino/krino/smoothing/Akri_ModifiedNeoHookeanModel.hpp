// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_MODIFIEDNEOHOOKEANMODEL_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_MODIFIEDNEOHOOKEANMODEL_HPP_

#include <limits>
#include <utility>

#include <Akri_KinematicUtils.hpp>
#include <MiniTensor.h>

namespace krino {

// Modified neo-Hookean energy density W(F, J) and its first Piola-Kirchhoff stress P = dW/dF, in 3D and 2D.
class ModifiedNeoHookeanModel
{
public:
  using Tensor = minitensor::Tensor<double, 3>;
  using Tensor2d = minitensor::Tensor<double, 2>;

  double energy(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const double Jm23 = kinematicUtils::pow_m23(J);
    const Tensor Cbar = Jm23 * C;

    const double WShape = minitensor::trace(Cbar) - 3.;

    const double Jn = 1.0 / J;
    const double WSize  = 1./4. * ((J - 1.0) * (J - 1.0) + (Jn - 1.0) * (Jn - 1.0));

    return WShape + WSize;
  }

  Tensor stress(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const Tensor FinvT = minitensor::transpose(minitensor::inverse(F));

    const double Jm23 = kinematicUtils::pow_m23(J);
    const Tensor Cbar = Jm23 * C;
    const double trCbar   = minitensor::trace(Cbar);

    const Tensor dWShape = 2. * Jm23 * F - 2./3. * trCbar * FinvT;

    const double Jn = 1.0 / J;
    const double J2  = J * J;
    const double Jn2 = 1.0 / J2;

    const Tensor dWSize = 1./2. * (J2 - J - Jn2 + Jn) * FinvT;

    return dWShape + dWSize;
  }

  double energy(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const double Jinv = 1./J;
    const Tensor2d Cbar = Jinv * C;

    const double WShape = minitensor::trace(Cbar) - 2.;

    const double Jn = 1.0 / J;
    const double WSize  = 1./4. * ((J - 1.0) * (J - 1.0) + (Jn - 1.0) * (Jn - 1.0));
    return WShape + WSize;
  }

  Tensor2d stress(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const Tensor2d FinvT = minitensor::transpose(minitensor::inverse(F));

    const double Jinv = 1./J;
    const Tensor2d Cbar = Jinv * C;
    const double trCbar = minitensor::trace(Cbar);

    const Tensor2d dWShape = 2.0 * Jinv * F - trCbar * FinvT;

    const double Jn = 1.0 / J;
    const double J2  = J * J;
    const double Jn2 = 1.0 / J2;

    const Tensor2d dWSize = 1./2. * (J2 - J - Jn2 + Jn) * FinvT;

    return dWShape + dWSize;
  }
};

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_MODIFIEDNEOHOOKEANMODEL_HPP_ */

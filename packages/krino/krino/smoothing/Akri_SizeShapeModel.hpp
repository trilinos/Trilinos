// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_SIZESHAPEMODEL_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_SIZESHAPEMODEL_HPP_

#include <limits>
#include <utility>

#include <Akri_KinematicUtils.hpp>
#include <MiniTensor.h>

namespace krino {

// Size-shape pseudo-energy density W(F, J) and its pseudo-stress P = dW/dF, in 3D and 2D.
class SizeShapeModel
{
public:
  using Tensor = minitensor::Tensor<double, 3>;
  using Tensor2d = minitensor::Tensor<double, 2>;

  double energy(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const double Jm23 = kinematicUtils::pow_m23(J);
    const Tensor Cbar = Jm23 * C;
    const double trCbar   = minitensor::trace(Cbar);

    const double etaShape = 1.0/3.0 * trCbar;
    const double etaSize = kinematicUtils::pow_p23(0.5*(J + 1./J));
    return etaShape*etaSize - 1.;
  }

  Tensor stress(const Tensor & F, const double J) const
  {
    const Tensor C = minitensor::transpose(F) * F;
    const Tensor FinvT = minitensor::transpose(minitensor::inverse(F));

    const double Jm23 = kinematicUtils::pow_m23(J);
    const Tensor Cbar = Jm23 * C;
    const double trCbar   = minitensor::trace(Cbar);

    const double etaShape = 1.0/3.0 * trCbar;
    const double etaSize = kinematicUtils::pow_p23(0.5*(J + 1./J));

    const Tensor dEtaShape = 2./3. * Jm23 * F - 2./9.*trCbar * FinvT;
    const Tensor dEtaSize = 1./3. * (J - 1./J) / std::cbrt(0.5*(J + 1./J)) * FinvT;

    return etaSize*dEtaShape + etaShape*dEtaSize;
  }

  double energy(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const double Jinv = 1./J;
    const Tensor2d Cbar = Jinv * C;
    const double trCbar   = minitensor::trace(Cbar);

    const double etaShape = 0.5 * trCbar;
    const double etaSize = 0.5*(J + Jinv);
    return etaShape*etaSize - 1.;
  }

  Tensor2d stress(const Tensor2d & F, const double J) const
  {
    const Tensor2d C = minitensor::transpose(F) * F;
    const Tensor2d FinvT = minitensor::transpose(minitensor::inverse(F));
    const double Jinv = 1./J;
    const Tensor2d Cbar = Jinv * C;
    const double trCbar = minitensor::trace(Cbar);

    const double etaShape = 0.5 * trCbar;
    const double etaSize = 0.5*(J + Jinv);

    const Tensor2d dEtaShape = Jinv * F - 0.5*trCbar * FinvT;
    const Tensor2d dEtaSize = 0.5 * (J - Jinv) * FinvT;

    return etaSize*dEtaShape + etaShape*dEtaSize;
  }
};

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_SIZESHAPEMODEL_HPP_ */

/*
 * Akri_PatchInterpolator.cpp
 *
 *  Created on: Feb 27, 2023
 *      Author: drnoble
 */
#include <Akri_PatchInterpolator.hpp>
#include <stk_math/StkVector.hpp>
#include <iostream>


namespace krino {

CubicSplineInterpolator::CubicSplineInterpolator(const std::array<stk::math::Vector3d,2> & sideNodeCoords, const std::array<stk::math::Vector3d,2> & sideNodeNormals)
{
  stk::math::Vector3d tangent = (sideNodeCoords[1]-sideNodeCoords[0]);
  const double dx = tangent.unitize();
  myNormal = crossZ(tangent);

  const double nDott_0 = -Dot(sideNodeNormals[0], tangent);
  const double y1_0 = nDott_0 / std::sqrt(1.-nDott_0*nDott_0);
  const double nDott_1 = -Dot(sideNodeNormals[1], tangent);
  const double y1_1 = nDott_1 / std::sqrt(1.-nDott_1*nDott_1);
  myY2H2Over6[0] = dx * (-2.*y1_0 - y1_1) / 3.;
  myY2H2Over6[1] = dx * (y1_0 + 2.*y1_1) / 3.;
}

double CubicSplineInterpolator::compute_spline_value(const double paramX) const
{
  const double a = 1.-paramX;
  const double b = paramX;

  return (a*a*a-a) * myY2H2Over6[0] + (b*b*b-b) * myY2H2Over6[1];
}

stk::math::Vector3d CubicSplineInterpolator::evaluate(const std::array<stk::math::Vector3d,2> & sideNodeCoords, const double paramX) const
{
  if (paramX == 0.)
    return sideNodeCoords[0];
  if (paramX == 1.)
    return sideNodeCoords[1];

  return (1.-paramX) * sideNodeCoords[0] + paramX * sideNodeCoords[1] + compute_spline_value(paramX) * myNormal;
}

}



/*
 * Akri_SimplexCalcs.cpp
 *
 *  Created on: Sep 11, 2025
 *      Author: drnoble
 */
#include <Akri_SimplexGradient.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

void calculate_triangle2d_gradient_and_area(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d, stk::math::Vector3d & grad, double & area)
{
  const double d10 = d[1] - d[0];
  const double d20 = d[2] - d[0];
  const stk::math::Vector3d x10 = x[1] - x[0];
  const stk::math::Vector3d x20 = x[2] - x[0];

  const double detJ = (x10[0]*x20[1]-x20[0]*x10[1]);
  area = detJ/2.;
  grad = (1./detJ)*(d10*stk::math::Vector3d(x20[1],-x20[0],0.0) + d20*stk::math::Vector3d(-x10[1],x10[0],0.0));
}

stk::math::Vector3d calculate_triangle2d_gradient(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d)
{
  stk::math::Vector3d grad;
  double area;
  calculate_triangle2d_gradient_and_area(x,d,grad,area);
  return grad;
}

double calculate_triangle2d_gradient_magnitude(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d)
{
  return calculate_triangle2d_gradient(x,d).length();
}

void calculate_tetrahedron_gradient_and_volume(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d, stk::math::Vector3d & grad, double & vol)
{
  const double d10 = d[1] - d[0];
  const double d20 = d[2] - d[0];
  const double d30 = d[3] - d[0];
  const stk::math::Vector3d x10 = x[1] - x[0];
  const stk::math::Vector3d x20 = x[2] - x[0];
  const stk::math::Vector3d x30 = x[3] - x[0];

  const stk::math::Vector3d x10_x_x20 = Cross(x10,x20);
  const stk::math::Vector3d x20_x_x30 = Cross(x20,x30);
  const stk::math::Vector3d x30_x_x10 = Cross(x30,x10);

  const double detJ = Dot(x30,x10_x_x20);
  vol = detJ/6.;
  grad = (1./detJ)*(d10*x20_x_x30 + d20*x30_x_x10 + d30*x10_x_x20);
}

stk::math::Vector3d calculate_tetrahedron_gradient(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d)
{
  stk::math::Vector3d grad;
  double vol;
  calculate_tetrahedron_gradient_and_volume(x,d,grad,vol);
  return grad;
}

double calculate_tetrahedron_gradient_magnitude(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d)
{
  return calculate_tetrahedron_gradient(x,d).length();
}

}



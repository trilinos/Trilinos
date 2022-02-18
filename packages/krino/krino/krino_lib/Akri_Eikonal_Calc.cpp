// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Eikonal_Calc.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Vec.hpp>
#include "Akri_MeshHelpers.hpp"

namespace krino {

double calculate_gradient_magnitude_triangle(const std::array<Vector3d,3> & x, const std::array<double,3> & d)
{
  const double d10 = d[1] - d[0];
  const double d20 = d[2] - d[0];
  const Vector3d x10 = x[1] - x[0];
  const Vector3d x20 = x[2] - x[0];

  const double detJ = (x10[0]*x20[1]-x20[0]*x10[1]);
  const Vector3d grad = d10*Vector3d(x20[1],-x20[0],0.0) + d20*Vector3d(-x10[1],x10[0],0.0);
  return grad.length()/detJ;
}

double calculate_gradient_magnitude_tetrahedron(const std::array<Vector3d,4> & x, const std::array<double,4> & d)
{
  const double d10 = d[1] - d[0];
  const double d20 = d[2] - d[0];
  const double d30 = d[3] - d[0];
  const Vector3d x10 = x[1] - x[0];
  const Vector3d x20 = x[2] - x[0];
  const Vector3d x30 = x[3] - x[0];

  const Vector3d x10_x_x20 = Cross(x10,x20);
  const Vector3d x20_x_x30 = Cross(x20,x30);
  const Vector3d x30_x_x10 = Cross(x30,x10);

  const double detJ = Dot(x30,x10_x_x20);
  const Vector3d grad = d10*x20_x_x30 + d20*x30_x_x10 + d30*x10_x_x20;
  return grad.length()/detJ;
}

double calculate_gradient_magnitude(const int npe,
    const stk::mesh::Entity * elem_nodes,
    const FieldRef dRef,
    const std::function<const Vector3d &(stk::mesh::Entity)> & get_coordinates)
{
  double mag_grad = 1.0;

  if (3 == npe)
  {
    const std::array<Vector3d,3> x{get_coordinates(elem_nodes[0]), get_coordinates(elem_nodes[1]), get_coordinates(elem_nodes[2])};
    const std::array<double,3> d{*field_data<double>(dRef, elem_nodes[0]), *field_data<double>(dRef, elem_nodes[1]), *field_data<double>(dRef, elem_nodes[2])};
    mag_grad = calculate_gradient_magnitude_triangle(x,d);
  }
  else
  {
    ThrowAssert(4 == npe);

    const std::array<Vector3d,4> x{get_coordinates(elem_nodes[0]), get_coordinates(elem_nodes[1]), get_coordinates(elem_nodes[2]), get_coordinates(elem_nodes[3])};
    const std::array<double,4> d{*field_data<double>(dRef, elem_nodes[0]), *field_data<double>(dRef, elem_nodes[1]), *field_data<double>(dRef, elem_nodes[2]), *field_data<double>(dRef, elem_nodes[3])};
    mag_grad = calculate_gradient_magnitude_tetrahedron(x,d);
  }

  return mag_grad;
}

std::array<int,3> get_oriented_nodes_triangle(const int nodeToUpdate)
{
  const std::array<int,3> lnn{(nodeToUpdate + 1) % 3, (nodeToUpdate + 2) % 3, nodeToUpdate};
  return lnn;
}

std::array<int,4> get_oriented_nodes_tetrahedron(const int nodeToUpdate)
{
  static const std::array<std::array<int,4>,4> nodePermute{{ {{1,3,2,0}}, {{0,2,3,1}}, {{0,3,1,2}}, {{0,1,2,3}} }};

  const std::array<int,4> lnn {nodePermute[nodeToUpdate][0], nodePermute[nodeToUpdate][1], nodePermute[nodeToUpdate][2], nodePermute[nodeToUpdate][3]};
  return lnn;
}

double eikonal_solve_triangle(const std::array<Vector3d,3> & x, const std::array<double,3> & d, const int sign, const int dim, const double speed)
{
  double dist_2 = d[2];

  if (sign*(dist_2-d[0]) < 0 || sign*(dist_2-d[1]) < 0)
  {
    return dist_2;
  }

  const double sqr_speed = speed*speed;
  const double d10 = d[1] - d[0];
  const Vector3d x10 = x[1] - x[0];
  const Vector3d x20 = x[2] - x[0];
  const double h10 = x10.length();
  const double h20 = x20.length();

  double detJ = 0;
  if (2 == dim)
  {
    detJ = (x10[0]*x20[1]-x20[0]*x10[1]);
  }
  else // (3 == dim)
  {
    detJ = Cross(x10,x20).length();
  }
  ThrowAssert(detJ > 0.0);

  const double a = h10*h10;
  const double b = -2.0 * sign*d10 * Dot(x10,x20);
  const double c = d10*d10 * h20*h20 - detJ*detJ/sqr_speed;

  bool elem_is_defining = false;

  const double det = b*b-4.0*a*c;
  if (det > 0.0)
  {
    // solve quadratic equation, roots are q/a and c/q
    const int sign_b = ( b < 0.0 ) ? -1 : 1;
    const double q = -0.5*(b + sign_b*std::sqrt(det));

    const double d20 = sign*std::max(c/q,q/a);

    const bool causal = (sign*d20 > 0.0 && sign*(d20-d10) > 0.0);

    if (causal)
    {
      const double loc = (Dot(x10,x20) - sqr_speed*d10*d20) / (h10*h10 - sqr_speed*d10*d10);
      elem_is_defining = (loc > 0.0 && loc < 1.0);

      if (elem_is_defining)
      {
        dist_2 = sign * std::min(sign*dist_2,sign*(d[0] + d20));
      }
    }
  }

  if (!elem_is_defining)
  {
    const double h21 = (x[2] - x[1]).length();
    dist_2 = sign * std::min(sign*dist_2,std::min(sign*d[0]+h20/speed,sign*d[1]+h21/speed));
    // Enforce causality - This is to catch the corner case (literally) where the characteristic is marching along the edges of the element
    dist_2 = sign * std::max(sign*dist_2,std::max(sign*d[0],sign*d[1]));
  }

  ThrowAssert(sign*(dist_2-d[0])>=0 && sign*(dist_2-d[1])>=0);

  return dist_2;
}

double eikonal_solve_tetrahedron(const std::array<Vector3d,4> & x, const std::array<double,4> & d, const int sign, const double speed)
{
  double dist_3 = d[3];
  if (sign*(dist_3-d[0]) < 0 || sign*(dist_3-d[1]) < 0 || sign*(dist_3-d[2]) < 0)
  {
    return dist_3;
  }

  const double sqr_speed = speed*speed;
  const double d10 = d[1]-d[0];
  const double d20 = d[2]-d[0];

  const Vector3d x10 = x[1] - x[0];
  const Vector3d x20 = x[2] - x[0];
  const Vector3d x30 = x[3] - x[0];

  const Vector3d x10_x_x20 = Cross(x10,x20);
  const Vector3d x20_x_x30 = Cross(x20,x30);
  const Vector3d x30_x_x10 = Cross(x30,x10);

  const double detJ = Dot(x30,x10_x_x20);
  ThrowAssert(detJ > 0);

  const Vector3d contrib12 = sign * (d10*x20_x_x30 + d20*x30_x_x10);

  const double a = x10_x_x20.length_squared();
  const double b = 2.0 * Dot(x10_x_x20,contrib12);
  const double c = contrib12.length_squared() - detJ*detJ/sqr_speed;

  bool elem_is_defining = false;

  const double det = b*b-4.0*a*c;
  if (det > 0.0)
  {
    // solve quadratic equation, roots are q/a and c/q
    const int sign_b = ( b < 0.0 ) ? -1 : 1;
    const double q = -0.5*(b + sign_b*std::sqrt(det));

    const double d30 = sign*std::max(c/q,q/a);

    const bool causal = (sign*d30 > 0.0 && sign*(d30-d10) > 0.0 && sign*(d30-d20) > 0.0);

    if (causal)
    {
      // Solve 2x2 system for parametric coords of intersection of gradient and 0-1-2
      // A1 r + B1 s == C1
      // A2 r + B2 s == C2
      const double A1 = x10.length_squared() - sqr_speed*d10*d10;
      const double B1 = Dot(x10,x20) - sqr_speed*d10*d20;
      const double C1 = Dot(x10,x30) - sqr_speed*d10*d30;
      const double A2 = B1;
      const double B2 = x20.length_squared() - sqr_speed*d20*d20;
      const double C2 = Dot(x20,x30) - sqr_speed*d20*d30;
      const double denom = A2*B1 - A1*B2;
      const double loc_r = (C2*B1 - C1*B2)/denom;
      const double loc_s = (A2*C1 - A1*C2)/denom;
      const double loc_t = 1.0 - loc_r - loc_s;
      elem_is_defining = (loc_r > 0.0 && loc_s > 0.0 && loc_t > 0.0);

      if (elem_is_defining)
      {
        dist_3 = sign * std::min(sign*dist_3,sign*(d[0] + d30));
      }
    }
  }

  if (!elem_is_defining)
  {
    static const std::array<std::array<int,3>,3> sides{{ {{0,1,3}}, {{1,2,3}}, {{2,0,3}} }};
    static constexpr int dim = 3;
    for (const auto & sideLNN : sides)
    {
      const std::array<Vector3d,3> xSide{x[sideLNN[0]],x[sideLNN[1]],x[sideLNN[2]]};
      const std::array<double,3> dSide{d[sideLNN[0]],d[sideLNN[1]],d[sideLNN[2]]};

      dist_3 = sign * std::min(sign*dist_3, sign*eikonal_solve_triangle(xSide,dSide,sign,dim,speed));
    }
    // Enforce causality - This is to catch the corner case (literally) where the characteristic is marching along the edges of the element
    // This is not completely covered by the 2d check because 1 face might still predict a value that is less than the neighbors.
    dist_3 = sign * std::max(sign*dist_3,std::max(sign*d[2],std::max(sign*d[1],sign*d[0])));
  }

  ThrowAssert(sign*(dist_3-d[0])>=0 && sign*(dist_3-d[1])>=0 && sign*(dist_3-d[2])>=0);

  return dist_3;
}

}

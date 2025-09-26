// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "Sacado.hpp"

#include <Akri_TriangleWithSensitivities.hpp>
#include <Akri_Triangle.hpp>

#include <stk_math/StkVector.hpp>

namespace krino {
namespace TriangleWithSens {

using TriSens = Sacado::Fad::SFad<double,9>;
using Vector3TriSens = stk::math::Vec<TriSens,3>;

Vector3TriSens Vector3TriSens_init(const stk::math::Vector3d & v, const int i)
{
  return Vector3TriSens(TriSens(9,3*i,v[0]), TriSens(9,3*i+1,v[1]), TriSens(9,3*i+2,v[2]));
}

Vector3TriSens Vector3TriSens_edge_init(const int i0, const stk::math::Vector3d & v0, const int i1, const stk::math::Vector3d & v1)
{
  Vector3TriSens vDx(TriSens(v1[0]-v0[0]), TriSens(v1[1]-v0[1]), TriSens(v1[2]-v0[2]));
  vDx[0].fastAccessDx(3*i0+0) = -1.;
  vDx[1].fastAccessDx(3*i0+1) = -1.;
  vDx[2].fastAccessDx(3*i0+2) = -1.;
  vDx[0].fastAccessDx(3*i1+0) = 1.;
  vDx[1].fastAccessDx(3*i1+1) = 1.;
  vDx[2].fastAccessDx(3*i1+2) = 1.;
  return vDx;
}

double area_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dArea)
{
  if (dArea == nullptr)
    return CalcTriangle3<double>::area(v0, v1, v2);

  const Vector3TriSens edgeDx0 = Vector3TriSens_edge_init(0,v0, 1,v1);
  const Vector3TriSens edgeDx1 = Vector3TriSens_edge_init(0,v0, 2,v2);

  const TriSens areaDx = CalcTriangle3<TriSens>::area(edgeDx0, edgeDx1);

  for (unsigned j=0; j<9; ++j)
    dArea[j] = areaDx.fastAccessDx(j);

  return areaDx.val();
}

stk::math::Vector3d normal_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dNormal)
{
  if (dNormal == nullptr)
    return CalcTriangle3<double>::normal(v0, v1, v2);

  const Vector3TriSens edgeDx0 = Vector3TriSens_edge_init(0,v0, 1,v1);
  const Vector3TriSens edgeDx1 = Vector3TriSens_edge_init(0,v0, 2,v2);

  const Vector3TriSens normalDx = CalcTriangle3<TriSens>::normal(edgeDx0, edgeDx1);

  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<9; ++j)
      dNormal[j*3+i] = normalDx[i].fastAccessDx(j);

  return stk::math::Vector3d(normalDx[0].val(), normalDx[1].val(), normalDx[2].val());
}

stk::math::Vector3d area_vector_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dAreaVector)
{
  if (dAreaVector == nullptr)
    return CalcTriangle3<double>::area_vector(v0, v1, v2);

  const Vector3TriSens edgeDx0 = Vector3TriSens_edge_init(0,v0, 1,v1);
  const Vector3TriSens edgeDx1 = Vector3TriSens_edge_init(0,v0, 2,v2);

  const Vector3TriSens areaVectorDx = CalcTriangle3<TriSens>::area_vector(edgeDx0, edgeDx1);

  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<9; ++j)
      dAreaVector[j*3+i] = areaVectorDx[i].fastAccessDx(j);

  return stk::math::Vector3d(areaVectorDx[0].val(), areaVectorDx[1].val(), areaVectorDx[2].val());
}

} // namespace TriangleWithSens
} // namespace krino

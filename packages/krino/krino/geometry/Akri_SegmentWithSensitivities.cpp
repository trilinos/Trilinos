// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "Sacado.hpp"

#include <Akri_SegmentWithSensitivities.hpp>
#include <Akri_Segment.hpp>

#include <stk_math/StkVector.hpp>

namespace krino {
namespace SegmentWithSens {

static constexpr unsigned numSens = 6;
using SegSens = Sacado::Fad::SFad<double,numSens>;
using Vector3SegSens = stk::math::Vec<SegSens,3>;


Vector3SegSens Vector3SegSens_edge_init(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1)
{
  Vector3SegSens vDx(SegSens(v1[0]-v0[0]), SegSens(v1[1]-v0[1]), SegSens(v1[2]-v0[2]));
  vDx[0].fastAccessDx(0) = -1.;
  vDx[1].fastAccessDx(1) = -1.;
  vDx[2].fastAccessDx(2) = -1.;
  vDx[0].fastAccessDx(3) = 1.;
  vDx[1].fastAccessDx(4) = 1.;
  vDx[2].fastAccessDx(5) = 1.;
  return vDx;
}

double length_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dLength)
{
  if (dLength == nullptr)
    return CalcSegment3<double>::length(v0, v1);

  const Vector3SegSens edgeDx = Vector3SegSens_edge_init(v0, v1);

  const SegSens lengthDx = edgeDx.length();

  for (unsigned j=0; j<numSens; ++j)
    dLength[j] = lengthDx.fastAccessDx(j);

  return lengthDx.val();
}

stk::math::Vector3d normal2d_and_optional_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dNormal)
{
  // This is the normal or the plane that passes through the 2 points and is perpendicular to the x-y plane
  if (dNormal == nullptr)
    return (crossZ(v1-v0)).unit_vector();

  const Vector3SegSens edgeDx = Vector3SegSens_edge_init(v0, v1);

  const Vector3SegSens normalDx = (crossZ(edgeDx)).unit_vector();

  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<numSens; ++j)
      dNormal[j*3+i] = normalDx[i].fastAccessDx(j);

  return stk::math::Vector3d(normalDx[0].val(), normalDx[1].val(), 0.);
}

} // namespace SegmentWithSens
} // namespace krino

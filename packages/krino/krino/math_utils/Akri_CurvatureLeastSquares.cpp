// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <vector>

#include <stk_util/util/SortAndUnique.hpp>
#include <Akri_CramersRuleSolver.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

static std::vector<int> get_unique_halo_nodes(const std::vector<std::array<int,2>> & haloSegments)
{
  std::vector<int> uniqueHaloNodes;
  for (auto && haloSegment : haloSegments)
  {
    uniqueHaloNodes.push_back(haloSegment[0]);
    uniqueHaloNodes.push_back(haloSegment[1]);
  }

  stk::util::sort_and_unique(uniqueHaloNodes);

  return uniqueHaloNodes;
}

void set_rotation_matrix_for_rotating_normal_to_zDir(std::array<std::array<double,3>,3> & m, const stk::math::Vector3d & normalDir)
{
  const stk::math::Vector3d normal = normalDir.unit_vector();
  static const stk::math::Vector3d zDir(0.,0.,1.);
  const double c = Dot(zDir, normal);
  stk::math::Vector3d v = Cross(normal, zDir);
  const double s = v.length();
  if (s > 0.) v *= (1./s);

  const double c1 = 1.-c;

  m[0][0] = c + v[0]*v[0]*c1;
  m[0][1] = v[0]*v[1]*c1 - v[2]*s;
  m[0][2] = v[0]*v[2]*c1 + v[1]*s;
  m[1][0] = v[1]*v[0]*c1 + v[2]*s;
  m[1][1] = c + v[1]*v[1]*(1.-c);
  m[1][2] = v[1]*v[2]*c1 - v[0]*s;
  m[2][0] = v[2]*v[0]*c1 - v[1]*s;
  m[2][1] = v[2]*v[1]*c1 + v[0]*s;
  m[2][2] = c + v[2]*v[2]*c1;
}

stk::math::Vector3d compute_patch_normal(const std::vector<stk::math::Vector3d> & haloNodeLocs, const std::vector<std::array<int,2>> & haloSegments)
{
  stk::math::Vector3d patchNormal = stk::math::Vector3d::ZERO;
  for (auto && haloSegment : haloSegments)
  {
    const stk::math::Vector3d & xc0 = haloNodeLocs[haloSegment[0]];
    const stk::math::Vector3d & xc1 = haloNodeLocs[haloSegment[1]];
    const stk::math::Vector3d wtNormal = Cross(xc0, xc1) / (xc0.length_squared()*xc1.length_squared());
    patchNormal += wtNormal;
  }

  return patchNormal.unit_vector();
}

static void fill_matrix_and_rhs_for_curvature_least_squares(const std::vector<stk::math::Vector3d> & rotatedUniqueHaloNodeLocs, std::array<std::array<double,3>,3> & A, std::array<double,3> & b)
{
  if (rotatedUniqueHaloNodeLocs.size() == 3)
  {
    for (int i=0; i<3; ++i)
    {
      A[i][0] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][0];
      A[i][1] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][1];
      A[i][2] = rotatedUniqueHaloNodeLocs[i][1]*rotatedUniqueHaloNodeLocs[i][1];
      b[i] = rotatedUniqueHaloNodeLocs[i][2];
    }
  }
  else
  {
    STK_ThrowRequireMsg(rotatedUniqueHaloNodeLocs.size() == 4, "Unexpected vector size in fill_matrix_and_rhs_for_curvature_least_squares.");
    std::array<std::array<double,3>,4> Apts;
    for (int i=0; i<4; ++i)
    {
      Apts[i][0] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][0];
      Apts[i][1] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][1];
      Apts[i][2] = rotatedUniqueHaloNodeLocs[i][1]*rotatedUniqueHaloNodeLocs[i][1];
    }

    for (int i=0; i<3; ++i)
    {
      b[i] = 0.;
      for (int k=0; k<4; ++k)
        b[i] += Apts[k][i] * rotatedUniqueHaloNodeLocs[k][2];

      for (int j=0; j<3; ++j)
      {
        A[i][j] = 0.;
        for (int k=0; k<4; ++k)
          A[i][j] += Apts[k][j]*Apts[k][i];
      }
    }
  }
}

static void fill_matrix_and_rhs_for_curvature_normal_least_squares(const std::vector<stk::math::Vector3d> & rotatedUniqueHaloNodeLocs, std::array<std::array<double,5>,5> & A, std::array<double,5> & b)
{
  if (rotatedUniqueHaloNodeLocs.size() == 5)
  {
    for (int i=0; i<5; ++i)
    {
      A[i][0] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][0];
      A[i][1] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][1];
      A[i][2] = rotatedUniqueHaloNodeLocs[i][1]*rotatedUniqueHaloNodeLocs[i][1];
      A[i][3] = rotatedUniqueHaloNodeLocs[i][0];
      A[i][4] = rotatedUniqueHaloNodeLocs[i][1];
      b[i] = rotatedUniqueHaloNodeLocs[i][2];
    }
  }
  else
  {
    std::vector<std::array<double,5>> Apts;
    Apts.resize(rotatedUniqueHaloNodeLocs.size());
    for (unsigned i=0; i<Apts.size(); ++i)
    {
      Apts[i][0] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][0];
      Apts[i][1] = rotatedUniqueHaloNodeLocs[i][0]*rotatedUniqueHaloNodeLocs[i][1];
      Apts[i][2] = rotatedUniqueHaloNodeLocs[i][1]*rotatedUniqueHaloNodeLocs[i][1];
      Apts[i][3] = rotatedUniqueHaloNodeLocs[i][0];
      Apts[i][4] = rotatedUniqueHaloNodeLocs[i][1];
    }

    for (unsigned i=0; i<5; ++i)
    {
      b[i] = 0.;
      for (unsigned k=0; k<Apts.size(); ++k)
        b[i] += Apts[k][i] * rotatedUniqueHaloNodeLocs[k][2];

      for (unsigned j=0; j<5; ++j)
      {
        A[i][j] = 0.;
        for (unsigned k=0; k<Apts.size(); ++k)
          A[i][j] += Apts[k][j]*Apts[k][i];
      }
    }
  }
}

static stk::math::Vector3d compute_least_squares_curvature_times_normal(const std::vector<stk::math::Vector3d> & rotatedUniqueHaloNodeLocs)
{
  if (rotatedUniqueHaloNodeLocs.size() == 3 || rotatedUniqueHaloNodeLocs.size() == 4)
  {
    std::array<std::array<double,3>,3> A;
    std::array<double,3> b;
    fill_matrix_and_rhs_for_curvature_least_squares(rotatedUniqueHaloNodeLocs, A, b);

    const std::array<double,3> soln = CramersRuleSolver::solve3x3(A,b);

    return stk::math::Vector3d(0., 0., -2.*soln[0]-2.*soln[2]);
  }
  else if (rotatedUniqueHaloNodeLocs.size() >= 5)
  {
    std::array<std::array<double,5>,5> A;
    std::array<double,5> b;
    fill_matrix_and_rhs_for_curvature_normal_least_squares(rotatedUniqueHaloNodeLocs, A, b);

    const std::array<double,5> soln = CramersRuleSolver::solve5x5(A,b);

    stk::math::Vector3d normal(-soln[3],-soln[4],1.);
    const double mag = normal.unitize();

    const double curvature =
        ((normal[0]*normal[0] - 1.) * 2.*soln[0] +
          normal[0]*normal[1] * 2.*soln[1] +
         (normal[1]*normal[1] - 1.) * 2.*soln[2]) / mag;

    return curvature*normal;
  }

  return stk::math::Vector3d::ZERO;
}

static stk::math::Vector3d compute_least_squares_normal(const std::vector<stk::math::Vector3d> & rotatedUniqueHaloNodeLocs)
{
  STK_ThrowRequire(rotatedUniqueHaloNodeLocs.size() >= 5);

  std::array<std::array<double,5>,5> A;
  std::array<double,5> b;
  fill_matrix_and_rhs_for_curvature_normal_least_squares(rotatedUniqueHaloNodeLocs, A, b);

  const std::array<double,5> soln = CramersRuleSolver::solve5x5(A,b);

  stk::math::Vector3d normal(-soln[3],-soln[4],1.);
  normal.unitize();

  return normal;
}

stk::math::Vector3d rotate_3d_vector(const std::array<std::array<double,3>,3> & m, const stk::math::Vector3d & v)
{
  return stk::math::Vector3d(
    (m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2]),
    (m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2]),
    (m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]));
}

stk::math::Vector3d reverse_rotate_3d_vector(const std::array<std::array<double,3>,3> & m, const stk::math::Vector3d & v)
{
  return stk::math::Vector3d(
    (m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2]),
    (m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2]),
    (m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2]));
}

static std::vector<stk::math::Vector3d> get_rotated_neighbor_node_locations(const std::vector<stk::math::Vector3d> & neighborNodeLocs, const std::array<std::array<double,3>,3> & m)
{
  std::vector<stk::math::Vector3d> rotatedUniqueHaloNodeLocs;
  rotatedUniqueHaloNodeLocs.reserve(neighborNodeLocs.size());
  for (auto && loc : neighborNodeLocs)
    rotatedUniqueHaloNodeLocs.push_back(rotate_3d_vector(m, loc));
  return rotatedUniqueHaloNodeLocs;
}

static std::vector<stk::math::Vector3d> get_rotated_unique_halo_node_locations(const std::vector<stk::math::Vector3d> & haloNodeLocs, std::vector<int> uniqueHaloNodes, const std::array<std::array<double,3>,3> & m)
{
  std::vector<stk::math::Vector3d> rotatedUniqueHaloNodeLocs;
  rotatedUniqueHaloNodeLocs.reserve(uniqueHaloNodes.size());
  for (int haloNode : uniqueHaloNodes)
    rotatedUniqueHaloNodeLocs.push_back(rotate_3d_vector(m, haloNodeLocs[haloNode]));
  return rotatedUniqueHaloNodeLocs;
}

stk::math::Vector3d compute_least_squares_curvature_times_normal(const std::vector<stk::math::Vector3d> & haloNodeLocs, const std::vector<std::array<int,2>> & haloSegments)
{
  if (haloSegments.size() < 3)
    return stk::math::Vector3d::ZERO;

  const stk::math::Vector3d patchNormal = compute_patch_normal(haloNodeLocs, haloSegments);

  std::vector<int> uniqueHaloNodes = get_unique_halo_nodes(haloSegments);

  std::array<std::array<double,3>,3> m;
  set_rotation_matrix_for_rotating_normal_to_zDir(m, patchNormal);

  const std::vector<stk::math::Vector3d> rotatedUniqueHaloNodeLocs = get_rotated_unique_halo_node_locations(haloNodeLocs, uniqueHaloNodes, m);

  const stk::math::Vector3d rotatedCurvatureNormal = compute_least_squares_curvature_times_normal(rotatedUniqueHaloNodeLocs);

  return reverse_rotate_3d_vector(m, rotatedCurvatureNormal);
}

stk::math::Vector3d compute_least_squares_curvature_times_normal(const stk::math::Vector3d & approximateNormal, const std::vector<stk::math::Vector3d> & neighborNodeLocs)
{
  if (neighborNodeLocs.size() < 3)
    return stk::math::Vector3d::ZERO;

  std::array<std::array<double,3>,3> m;
  set_rotation_matrix_for_rotating_normal_to_zDir(m, approximateNormal);

  const std::vector<stk::math::Vector3d> rotatedNbrNodeLocs = get_rotated_neighbor_node_locations(neighborNodeLocs, m);

  const stk::math::Vector3d rotatedCurvatureNormal = compute_least_squares_curvature_times_normal(rotatedNbrNodeLocs);

  return reverse_rotate_3d_vector(m, rotatedCurvatureNormal);
}

stk::math::Vector3d compute_least_squares_normal(const stk::math::Vector3d & approximateNormal, const std::vector<stk::math::Vector3d> & neighborNodeLocs)
{
  if (neighborNodeLocs.size() < 5)
    return approximateNormal;

  std::array<std::array<double,3>,3> m;
  set_rotation_matrix_for_rotating_normal_to_zDir(m, approximateNormal);

  const std::vector<stk::math::Vector3d> rotatedNbrNodeLocs = get_rotated_neighbor_node_locations(neighborNodeLocs, m);

  const stk::math::Vector3d rotatedCurvatureNormal = compute_least_squares_normal(rotatedNbrNodeLocs);

  return reverse_rotate_3d_vector(m, rotatedCurvatureNormal).unit_vector();
}

}



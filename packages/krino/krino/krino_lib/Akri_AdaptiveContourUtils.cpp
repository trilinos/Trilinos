/*
 * Akri_AdaptiveContourUtils.cpp
 *
 *  Created on: Sep 30, 2025
 *      Author: drnoble
 */

#include <Akri_AdaptiveContourUtils.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_Sign.hpp>

namespace krino {

template <size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NNODES> & coords,
    const std::array<double,NNODES> & dist)
{
  for (size_t n1=0; n1<NNODES; ++n1)
  {
    for (size_t n2=n1; n2<NNODES; ++n2)
    {
      const double sqrLen = (coords[n1] - coords[n2]).length_squared();
      if (dist[n1]*dist[n1] <= 0.25*sqrLen || dist[n2]*dist[n2] <= 0.25*sqrLen)
        return true;
    }
  }
  return false;
}

template <size_t NPARENTNODES, size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NPARENTNODES> & coords,
    const std::array<double,NPARENTNODES> & dist,
    const std::array<int, NNODES> & subElemVertexIndices)
{
  for (size_t n1=0; n1<NNODES; ++n1)
  {
    const int pn1 = subElemVertexIndices[n1];
    for (size_t n2=n1; n2<NNODES; ++n2)
    {
      const int pn2 = subElemVertexIndices[n2];
      const double sqrLen = (coords[pn1] - coords[pn2]).length_squared();
      if (dist[pn1]*dist[pn1] <= 0.25*sqrLen || dist[pn2]*dist[pn2] <= 0.25*sqrLen)
        return true;
    }
  }
  return false;
}

double clip_midedge_distance(const double d0, const double d1, const double d2)
{
  const double d25 = 0.75*d0+0.25*d1;
  const double d75 = 0.25*d0+0.75*d1;
  if (d0 < d1)
    return std::max(std::min(d2, d75), d25);
  return std::max(std::min(d2, d25), d75);
}

bool is_edge_converged(const double d0, const double d1, const double d2, const double nonlinearDistTol)
{
  return (std::abs(d2 - 0.5*(d0+d1)) < nonlinearDistTol);
}

int compute_node_sign(const double dist)
{
  return (dist == 0) ? 0 : ((dist < 0) ? -1 : 1);
}

template <size_t NCOORDS, size_t NDIST>
stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NDIST> & distance,
  const unsigned i0, const unsigned i1, const unsigned i2)
{
  const double loc = find_quadratic_crossing(distance[i0], distance[i1], distance[i2]);
  return (1.-loc) * coords[i0] + loc * coords[i1];
}

// Explicit template instantiation
template bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,3> & coords, const std::array<double,3> & dist);
template bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,4> & coords, const std::array<double,4> & dist);
template bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,6> & coords, const std::array<double,6> & dist, const std::array<int, 3> & subElemVertexIndices);
template bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,10> & coords, const std::array<double,10> & dist, const std::array<int, 4> & subElemVertexIndices);
template stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,3> & coords, const std::array<double,6> & distance, const unsigned i0, const unsigned i1, const unsigned i2);
template stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,4> & coords, const std::array<double,10> & distance, const unsigned i0, const unsigned i1, const unsigned i2);

}



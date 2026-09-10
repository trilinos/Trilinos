// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURUTILS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURUTILS_HPP_
#include <stk_math/StkVector.hpp>

#include <array>

namespace krino {

template <size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NNODES> & coords,
    const std::array<double,NNODES> & dist);

template <size_t NPARENTNODES, size_t NNODES>
bool have_possibly_cut_edge(const std::array<stk::math::Vector3d,NPARENTNODES> & coords,
    const std::array<double,NPARENTNODES> & dist,
    const std::array<int, NNODES> & subElemVertexIndices);

double clip_midedge_distance(const double d0, const double d1, const double d2);

bool is_edge_converged(const double d0, const double d1, const double d2, const double nonlinearDistTol);

int compute_node_sign(const double dist);

template <size_t NNODES>
void snap_distance(std::array<double,NNODES> & dist, const double snapDistTol)
{
  for (size_t n=0; n<NNODES; ++n)
    dist[n] = (std::abs(dist[n]) < snapDistTol) ? 0.0 : dist[n];
}

template <size_t NNODES, size_t NSUBNODES, typename T>
std::array<T, NSUBNODES> subarray(const std::array<T, NNODES> & a, const std::array<int,NSUBNODES> & indices)
{
  std::array<T, NSUBNODES> subarray;
  for (size_t n=0; n<NSUBNODES; ++n)
    subarray[n] = a[indices[n]];
  return subarray;
}

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_ADAPTIVECONTOURUTILS_HPP_ */

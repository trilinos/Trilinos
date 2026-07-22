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

template <size_t NCOORDS, size_t NDIST>
stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NDIST> & distance,
  const unsigned i0, const unsigned i1, const unsigned i2);

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

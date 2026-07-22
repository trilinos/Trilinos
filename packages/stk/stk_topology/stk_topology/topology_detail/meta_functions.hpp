// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STKTOPOLOGY_DETAIL_META_FUNCTION_HPP
#define STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

#include "stk_topology/topology_decl.hpp"
#include "stk_util/stk_config.h"
#include <type_traits>

namespace stk::topology_detail {

template <typename Topology, unsigned EdgeOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr unsigned num_edge_nodes_() {
  return (Topology::edge_node_ordinals_offsets[EdgeOrdinal+1] - Topology::edge_node_ordinals_offsets[EdgeOrdinal]);
}

template <typename Topology, unsigned FaceOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr unsigned num_face_nodes_() {
  return (Topology::face_node_ordinals_offsets[FaceOrdinal+1] - Topology::face_node_ordinals_offsets[FaceOrdinal]);
}

template <typename Topology, unsigned EdgeOrdinal, unsigned NodeOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr unsigned edge_node_ordinal_()
{
  return (Topology::edge_node_ordinals_vector[Topology::edge_node_ordinals_offsets[EdgeOrdinal] + NodeOrdinal]);
}

template <typename Topology, unsigned FaceOrdinal, unsigned NodeOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr unsigned face_node_ordinal_()
{
  return (Topology::face_node_ordinals_vector[Topology::face_node_ordinals_offsets[FaceOrdinal] + NodeOrdinal]);
}

template <typename Topology, unsigned PermutationOrdinal, unsigned NodeOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr unsigned permutation_node_ordinal_()
{
  return Topology::permutation_node_ordinals_vector[PermutationOrdinal*Topology::num_nodes + NodeOrdinal];
}


//------------------------------------------------------------------------------
template <typename Topology, unsigned SpatialDimension>
KOKKOS_INLINE_FUNCTION
constexpr bool defined_on_spatial_dimension_()
{
  static_assert(SpatialDimension < 4, "Invalid spatial dimension");
  return Topology::spatial_dimension_vector[SpatialDimension];
}

//------------------------------------------------------------------------------

template <typename Topology, unsigned EdgeOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr topology::topology_t edge_topology_()
{
  if constexpr (EdgeOrdinal < Topology::num_edges)
  {
    return Topology::edge_topology_vector[EdgeOrdinal]; 
  }
  return topology::INVALID_TOPOLOGY;
}

//------------------------------------------------------------------------------

template <typename Topology, unsigned FaceOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr topology::topology_t face_topology_()
{
  if constexpr (FaceOrdinal < Topology::num_faces)
  {
    return Topology::face_topology_vector[FaceOrdinal];
  }
  return topology::INVALID_TOPOLOGY;
}

//------------------------------------------------------------------------------

template <typename Topology, unsigned SideOrdinal>
KOKKOS_INLINE_FUNCTION
constexpr topology::rank_t side_rank_()
{
  if constexpr (SideOrdinal < Topology::num_faces) {
    return topology::FACE_RANK;
  } else {
    return topology::EDGE_RANK;
  }
  return Topology::side_rank;
}

template <typename Topology>
KOKKOS_INLINE_FUNCTION
constexpr unsigned num_side_ranks_() {
  if constexpr (Topology::has_mixed_rank_sides) {
    return 2u;
  }
  return (Topology::side_rank != topology::INVALID_RANK ? 1u : 0u);
}

template <typename Topology, typename SideRankOutputIterator>
KOKKOS_INLINE_FUNCTION
constexpr void side_ranks_( SideRankOutputIterator output_ranks )
{
  if constexpr (num_side_ranks_<Topology>() == 2) {
    *output_ranks = topology::FACE_RANK; ++output_ranks;
    *output_ranks = topology::EDGE_RANK;
  } else if constexpr (num_side_ranks_<Topology>() == 1) {
    *output_ranks = Topology::side_rank;
  }
}

//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned EdgeOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct edge_node_ordinals_impl_ {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(edge_node_ordinal_<Topology, EdgeOrdinal, CurrentNode>()),
           edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned EdgeOrdinal, unsigned NumNodes>
struct edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, NumNodes, NumNodes> {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor /*fillOutput*/) {
    return 0;
  }
};

template <typename Topology, unsigned EdgeOrdinal, typename OrdinalOutputFunctor>
KOKKOS_INLINE_FUNCTION
constexpr int edge_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  if constexpr (EdgeOrdinal < Topology::num_edges)
  {
    return edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, num_edge_nodes_<Topology, EdgeOrdinal>()>::execute(fillOutput);
  }
  return 0;
}


//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned FaceOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct face_node_ordinals_impl_ {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(face_node_ordinal_<Topology, FaceOrdinal, CurrentNode>()),
           face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned FaceOrdinal, unsigned NumNodes>
struct face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, NumNodes, NumNodes> {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor /*fillOutput*/) {
    return 0;
  }
};

template <typename Topology, unsigned FaceOrdinal, typename OrdinalOutputFunctor>
KOKKOS_INLINE_FUNCTION
constexpr int face_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  if constexpr (FaceOrdinal < Topology::num_faces)
  {
    return face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, num_face_nodes_<Topology, FaceOrdinal>()>::execute(fillOutput);
  }
  return 0;
}


//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned PermutationOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct permutation_node_ordinals_impl_ {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(permutation_node_ordinal_<Topology, PermutationOrdinal, CurrentNode>()),
           permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned PermutationOrdinal, unsigned NumNodes>
struct permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, NumNodes, NumNodes> {
  KOKKOS_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor /*fillOutput*/) {
    return 0;
  }
};

template <typename Topology, unsigned PermutationOrdinal, typename OrdinalOutputFunctor>
KOKKOS_INLINE_FUNCTION
constexpr int permutation_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  if constexpr (PermutationOrdinal < Topology::num_permutations)
  {
    return permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, Topology::num_nodes>::execute(fillOutput);
  }
  return 0;
}

//------------------------------------------------------------------------------
} //namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

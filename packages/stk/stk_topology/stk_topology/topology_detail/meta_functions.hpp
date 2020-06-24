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

#include <type_traits>

namespace stk { namespace topology_detail {

//------------------------------------------------------------------------------
template <typename Topology, unsigned EdgeOrdinal>
STK_INLINE_FUNCTION
constexpr unsigned num_edge_nodes_() {
  return (Topology::edge_node_ordinals_offsets[EdgeOrdinal+1] - Topology::edge_node_ordinals_offsets[EdgeOrdinal]);
}

template <typename Topology, unsigned FaceOrdinal>
STK_INLINE_FUNCTION
constexpr unsigned num_face_nodes_() {
  return (Topology::face_node_ordinals_offsets[FaceOrdinal+1] - Topology::face_node_ordinals_offsets[FaceOrdinal]);
}

template <typename Topology, unsigned EdgeOrdinal, unsigned NodeOrdinal>
STK_INLINE_FUNCTION
constexpr unsigned edge_node_ordinal_()
{
  return (Topology::edge_node_ordinals_vector[Topology::edge_node_ordinals_offsets[EdgeOrdinal] + NodeOrdinal]);
}

template <typename Topology, unsigned FaceOrdinal, unsigned NodeOrdinal>
STK_INLINE_FUNCTION
constexpr unsigned face_node_ordinal_()
{
  return (Topology::face_node_ordinals_vector[Topology::face_node_ordinals_offsets[FaceOrdinal] + NodeOrdinal]);
}

template <typename Topology, unsigned PermutationOrdinal, unsigned NodeOrdinal>
STK_INLINE_FUNCTION
constexpr unsigned permutation_node_ordinal_()
{
  return Topology::permutation_node_ordinals_vector[PermutationOrdinal*Topology::num_nodes + NodeOrdinal];
}


//------------------------------------------------------------------------------
template <typename Topology, unsigned SpatialDimension>
STK_INLINE_FUNCTION
constexpr bool defined_on_spatial_dimension_()
{
  static_assert(SpatialDimension < 4, "Invalid spatial dimension");
  return Topology::spatial_dimension_vector[SpatialDimension];
}

//------------------------------------------------------------------------------

template <typename Topology, unsigned EdgeOrdinal>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<EdgeOrdinal < Topology::num_edges, topology::topology_t>::type
edge_topology_()
{
  return Topology::edge_topology_vector[EdgeOrdinal];
}

template <typename Topology, unsigned EdgeOrdinal>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<EdgeOrdinal >= Topology::num_edges, topology::topology_t>::type
edge_topology_()
{
  return topology::INVALID_TOPOLOGY;
}

//------------------------------------------------------------------------------

template <typename Topology, unsigned FaceOrdinal>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<FaceOrdinal < Topology::num_faces, topology::topology_t>::type
face_topology_()
{
  return Topology::face_topology_vector[FaceOrdinal];
}

template <typename Topology, unsigned FaceOrdinal>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<FaceOrdinal >= Topology::num_faces, topology::topology_t>::type
face_topology_()
{
  return topology::INVALID_TOPOLOGY;
}

//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned EdgeOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct edge_node_ordinals_impl_ {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(edge_node_ordinal_<Topology, EdgeOrdinal, CurrentNode>()),
           edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned EdgeOrdinal, unsigned NumNodes>
struct edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, NumNodes, NumNodes> {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return 0;
  }
};

template <typename Topology, unsigned EdgeOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(EdgeOrdinal < Topology::num_edges), int>::type edge_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return edge_node_ordinals_impl_<Topology, OrdinalOutputFunctor, EdgeOrdinal, num_edge_nodes_<Topology, EdgeOrdinal>()>::execute(fillOutput);
}

template <typename Topology, unsigned EdgeOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(EdgeOrdinal >= Topology::num_edges), int>::type edge_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return 0;
}


//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned FaceOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct face_node_ordinals_impl_ {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(face_node_ordinal_<Topology, FaceOrdinal, CurrentNode>()),
           face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned FaceOrdinal, unsigned NumNodes>
struct face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, NumNodes, NumNodes> {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return 0;
  }
};

template <typename Topology, unsigned FaceOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(FaceOrdinal < Topology::num_faces), int>::type face_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return face_node_ordinals_impl_<Topology, OrdinalOutputFunctor, FaceOrdinal, num_face_nodes_<Topology, FaceOrdinal>()>::execute(fillOutput);
}

template <typename Topology, unsigned FaceOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(FaceOrdinal >= Topology::num_faces), int>::type face_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return 0;
}


//------------------------------------------------------------------------------
template <typename Topology, typename OrdinalOutputFunctor, unsigned PermutationOrdinal, unsigned NumNodes, unsigned CurrentNode = 0>
struct permutation_node_ordinals_impl_ {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return fillOutput(permutation_node_ordinal_<Topology, PermutationOrdinal, CurrentNode>()),
           permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, NumNodes, CurrentNode+1>::execute(fillOutput);
  }
};

template <typename Topology, typename OrdinalOutputFunctor, unsigned PermutationOrdinal, unsigned NumNodes>
struct permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, NumNodes, NumNodes> {
  STK_INLINE_FUNCTION
  constexpr static int execute(OrdinalOutputFunctor fillOutput) {
    return 0;
  }
};

template <typename Topology, unsigned PermutationOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(PermutationOrdinal < Topology::num_permutations), int>::type permutation_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return permutation_node_ordinals_impl_<Topology, OrdinalOutputFunctor, PermutationOrdinal, Topology::num_nodes>::execute(fillOutput);
}

template <typename Topology, unsigned PermutationOrdinal, typename OrdinalOutputFunctor>
STK_INLINE_FUNCTION
constexpr typename std::enable_if<(PermutationOrdinal >= Topology::num_permutations), int>::type permutation_node_ordinals_(OrdinalOutputFunctor fillOutput)
{
  return 0;
}


//------------------------------------------------------------------------------
}} //namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

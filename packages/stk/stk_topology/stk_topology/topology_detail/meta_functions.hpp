// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

#define STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(name,result)     \
  template <typename Topology>                                   \
  struct name##_                                       \
    : public boost::mpl::integral_c< result, Topology::name >    \
  {};

namespace stk { namespace topology_detail {

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(is_valid,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(is_shell,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(has_homogeneous_faces,bool)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(rank,topology::rank_t)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(side_rank,topology::rank_t)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(dimension,unsigned)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_nodes,unsigned)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_vertices,unsigned)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_edges,unsigned)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_faces,unsigned)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_permutations,unsigned)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(base,topology::topology_t)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(edge_topology,topology::topology_t)


template <typename Topology, unsigned SpatialDimension>
struct defined_on_spatial_dimension_
  : public boost::mpl::eval_if_c<  (SpatialDimension < 4)
      , boost::mpl::at_c< typename Topology::spatial_dimension_vector, SpatialDimension >
      , boost::mpl::identity< boost::mpl::false_>
    >::type
{};

template <typename Topology, unsigned EdgeOrdinal>
struct edge_node_ordinals_
  : public boost::mpl::eval_if_c<  (EdgeOrdinal < num_edges_<Topology>::value)
      , boost::mpl::at_c< typename Topology::edge_node_ordinals_vector, EdgeOrdinal >
      , boost::mpl::identity< boost::mpl::vector_c<unsigned> >
    >::type
{};

template <typename Topology, unsigned FaceOrdinal>
struct face_topology_
  : public boost::mpl::eval_if_c<  (FaceOrdinal < num_faces_<Topology>::value)
      , boost::mpl::at_c< typename Topology::face_topology_vector, FaceOrdinal >
      , boost::mpl::identity< boost::mpl::integral_c<topology::topology_t, topology::INVALID_TOPOLOGY> >
    >::type
{};

template <typename Topology, unsigned FaceOrdinal>
struct face_node_ordinals_
  : public boost::mpl::eval_if_c<  (FaceOrdinal < num_faces_<Topology>::value)
      , boost::mpl::at_c< typename Topology::face_node_ordinals_vector, FaceOrdinal >
      , boost::mpl::identity< boost::mpl::vector_c<unsigned> >
    >::type
{};

template <typename Topology, unsigned PermutationOrdinal>
struct permutation_node_ordinals_
  : public boost::mpl::eval_if_c<  (PermutationOrdinal < num_permutations_<Topology>::value)
      , boost::mpl::at_c< typename Topology::permutation_node_ordinals_vector, PermutationOrdinal >
      , boost::mpl::identity< boost::mpl::vector_c<unsigned> >
    >::type
{};

}} //namespace stk::topology_detail

#undef STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION

#endif //STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

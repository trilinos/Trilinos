// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

namespace stk { namespace topology_detail {

template <typename Topology>
struct is_valid_ : public boost::mpl::integral_c<bool, Topology::is_valid> {};

template <typename Topology>
struct is_shell_ : public boost::mpl::integral_c<bool, Topology::is_shell> {};

template <typename Topology>
struct has_homogeneous_faces_ : public boost::mpl::integral_c<bool, Topology::has_homogeneous_faces> {};

template <typename Topology>
struct rank_ : public boost::mpl::integral_c<topology::rank_t, Topology::rank> {};

template <typename Topology>
struct side_rank_ : public boost::mpl::integral_c<topology::rank_t, Topology::side_rank> {};

template <typename Topology>
struct dimension_ : public boost::mpl::integral_c<unsigned, Topology::dimension> {};

template <typename Topology>
struct num_nodes_ : public boost::mpl::integral_c<unsigned, Topology::num_nodes> {};

template <typename Topology>
struct num_vertices_ : public boost::mpl::integral_c<unsigned, Topology::num_vertices> {};

template <typename Topology>
struct num_edges_ : public boost::mpl::integral_c<unsigned, Topology::num_edges> {};

template <typename Topology>
struct num_faces_ : public boost::mpl::integral_c<unsigned, Topology::num_faces> {};

template <typename Topology>
struct num_permutations_ : public boost::mpl::integral_c<unsigned, Topology::num_permutations> {};

template <typename Topology>
struct base_ : public boost::mpl::integral_c<topology::topology_t, Topology::base> {};

template <typename Topology>
struct edge_topology_ : public boost::mpl::integral_c<topology::topology_t, Topology::edge_topology> {};

template <typename Topology, unsigned SpatialDimension>
struct defined_on_spatial_dimension_
  : public boost::mpl::eval_if_c< (SpatialDimension < 4),
                                  boost::mpl::at_c<typename Topology::spatial_dimension_vector, SpatialDimension>,
                                  boost::mpl::identity<boost::mpl::false_> >::type
{};

template <typename Topology, unsigned EdgeOrdinal>
struct edge_node_ordinals_
  : public boost::mpl::eval_if_c< (EdgeOrdinal < num_edges_<Topology>::value),
                                  boost::mpl::at_c<typename Topology::edge_node_ordinals_vector, EdgeOrdinal>,
                                  boost::mpl::identity<boost::mpl::vector_c<unsigned>> >::type
{};

template <typename Topology, unsigned FaceOrdinal>
struct face_topology_
  : public boost::mpl::eval_if_c< (FaceOrdinal < num_faces_<Topology>::value),
                                  boost::mpl::at_c<typename Topology::face_topology_vector, FaceOrdinal>,
                                  boost::mpl::identity<boost::mpl::integral_c<topology::topology_t, topology::INVALID_TOPOLOGY>> >::type
{};

template <typename Topology, unsigned FaceOrdinal>
struct face_node_ordinals_
  : public boost::mpl::eval_if_c< (FaceOrdinal < num_faces_<Topology>::value),
                                  boost::mpl::at_c<typename Topology::face_node_ordinals_vector, FaceOrdinal>,
                                  boost::mpl::identity<boost::mpl::vector_c<unsigned>> >::type
{};

template <typename Topology, unsigned PermutationOrdinal>
struct permutation_node_ordinals_
  : public boost::mpl::eval_if_c< (PermutationOrdinal < num_permutations_<Topology>::value),
                                  boost::mpl::at_c<typename Topology::permutation_node_ordinals_vector, PermutationOrdinal >,
                                  boost::mpl::identity<boost::mpl::vector_c<unsigned>> >::type
{};

}} //namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

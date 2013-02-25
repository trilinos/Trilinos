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
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(is_heterogenuous,bool)
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

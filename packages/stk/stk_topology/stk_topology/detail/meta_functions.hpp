#ifndef STKTOPOLOGY_DETAIL_META_FUNCTION_HPP
#define STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>

#define STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(name,result)     \
  template <typename Topology>                                   \
  struct name##_                                       \
    : public boost::mpl::integral_c< result, Topology::name >    \
  {};

namespace stk { namespace detail {

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(is_valid,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(is_shell,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(has_homogeneous_edges,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(has_homogeneous_faces,bool)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(has_homogeneous_sides,bool)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(rank,topology::rank_t)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(side_rank,topology::rank_t)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(dimension,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_nodes,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_vertices,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_edges,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_faces,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_sides,int)
STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(num_permutations,int)

STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION(base,topology::topology_t)


template <typename Topology, int SpatialDimension>
struct defined_on_spatial_dimension_
  : public boost::mpl::if_c<  (SpatialDimension < 4)
      , typename boost::mpl::at_c< typename Topology::spatial_dimension_vector, SpatialDimension >::type
      , boost::mpl::false_
    >::type
{};

template <typename Topology, int EdgeOrdinal>
struct edge_topology_
  : public boost::mpl::if_c<  (EdgeOrdinal < num_edges_<Topology>::value)
      , typename boost::mpl::at_c< typename Topology::edge_topology_vector, EdgeOrdinal >::type
      , boost::mpl::integral_c<topology::topology_t, topology::INVALID_TOPOLOGY>
    >::type
{};

template <typename Topology, int EdgeOrdinal>
struct edge_node_ordinals_
  : public boost::mpl::if_c<  (EdgeOrdinal < num_edges_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::edge_node_ordinals_vector, EdgeOrdinal >::type
      , boost::mpl::vector_c<int>
    >::type
{};

template <typename Topology, int FaceOrdinal>
struct face_topology_
  : public boost::mpl::if_c<  (FaceOrdinal < num_faces_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::face_topology_vector, FaceOrdinal >::type
      , boost::mpl::integral_c<topology::topology_t, topology::INVALID_TOPOLOGY>
    >::type
{};

template <typename Topology, int FaceOrdinal>
struct face_node_ordinals_
  : public boost::mpl::if_c<  (FaceOrdinal < num_faces_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::face_node_ordinals_vector, FaceOrdinal >::type
      , boost::mpl::vector_c<int>
    >::type
{};

template <typename Topology, int FaceOrdinal>
struct side_topology_
  : public boost::mpl::if_c<  (FaceOrdinal < num_sides_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::side_topology_vector, FaceOrdinal >::type
      , boost::mpl::integral_c<topology::topology_t, topology::INVALID_TOPOLOGY>
    >::type
{};

template <typename Topology, int FaceOrdinal>
struct side_node_ordinals_
  : public boost::mpl::if_c<  (FaceOrdinal < num_sides_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::side_node_ordinals_vector, FaceOrdinal >::type
      , boost::mpl::vector_c<int>
    >::type
{};

template <typename Topology, int PermutationOrdinal>
struct permutation_node_ordinals_
  : public boost::mpl::if_c<  (PermutationOrdinal < num_permutations_<Topology>::value)
      , typename  boost::mpl::at_c< typename Topology::permutation_node_ordinals_vector, PermutationOrdinal >::type
      , boost::mpl::vector_c<int>
    >::type
{};

}} //namespace stk::detail

#undef STKTOPOLOGY_DETAIL_SIMPLE_META_FUNCTION

#endif //STKTOPOLOGY_DETAIL_META_FUNCTION_HPP

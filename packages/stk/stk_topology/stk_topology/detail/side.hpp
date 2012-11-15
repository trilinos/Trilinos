#ifndef STKTOPOLOGY_DETAIL_SIDE_HPP
#define STKTOPOLOGY_DETAIL_SIDE_HPP

#include <stk_topology/topology.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>

namespace stk { namespace detail {
// mpl vector of side topologies
// NOTE: This assumes that the only topologies with a side rank of NODE_RANK is a LINE or PARTICLE
template <typename TopologyData>
struct side_topology_vector_helper
  : public
    boost::mpl::if_c< (topology::FACE_RANK == TopologyData::side_rank)
      , typename TopologyData::face_topology_vector
      , typename boost::mpl::if_c< (topology::EDGE_RANK == TopologyData::side_rank)
        , typename TopologyData::edge_topology_vector
        , typename boost::mpl::if_c< (topology::NODE_RANK == TopologyData::side_rank)
          , boost::mpl::vector_c<topology::topology_t, topology::NODE, topology::NODE>
          , boost::mpl::vector_c<topology::topology_t>
        >::type
      >::type
    >::type
{};

// NOTE: This assumes that the only topologies with a side rank of NODE_RANK is a LINE or PARTICLE
template <typename TopologyData>
struct side_node_ordinals_vector_helper
  : public
    boost::mpl::if_c< (topology::FACE_RANK == TopologyData::side_rank)
      , typename TopologyData::face_node_ordinals_vector
      , typename boost::mpl::if_c< (topology::EDGE_RANK == TopologyData::side_rank)
        , typename TopologyData::edge_node_ordinals_vector
        , typename boost::mpl::if_c< (topology::NODE_RANK == TopologyData::side_rank)
          , boost::mpl::vector< boost::mpl::vector_c<int, 0>, boost::mpl::vector_c<int, 1> >
          , boost::mpl::vector<>
        >::type
      >::type
    >::type
{};

template <typename TopologyData>
struct num_sides_helper
{
  typedef int value_type;
  static const int value =    (topology::FACE_RANK == TopologyData::side_rank) ? TopologyData::num_faces
                                : ((topology::EDGE_RANK == TopologyData::side_rank) ? TopologyData::num_edges
                                : ((topology::NODE_RANK == TopologyData::side_rank) ? TopologyData::num_vertices : 0) );
};

template <typename TopologyData>
struct has_homogeneous_sides_helper
{
  typedef bool value_type;
  static const bool value =    (topology::FACE_RANK == TopologyData::side_rank) ? TopologyData::has_homogeneous_faces
                            : ((topology::EDGE_RANK == TopologyData::side_rank) ? TopologyData::has_homogeneous_edges
                            :  (topology::NODE_RANK == TopologyData::side_rank) );
};


}} // namespace stk::detail

#endif //STKTOPOLOGY_DETAIL_SIDE_HPP

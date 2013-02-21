#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>


#include <boost/mpl/assert.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>

using namespace stk;
using namespace stk::topology_detail;
using namespace boost;

namespace {

template <typename TopologyData, typename PermutationVector, unsigned NumNodes, unsigned Permutation = 0, unsigned NumPermutations = mpl::size<PermutationVector>::value>
struct check_permutations
{
  static const bool value =    (NumNodes == mpl::size< typename mpl::at_c<PermutationVector,Permutation>::type>::value)
                            && check_permutations<TopologyData,PermutationVector,NumNodes,Permutation+1>::value;

  BOOST_MPL_ASSERT_MSG(   value
                        , PERMUTATION_NODE_ORDINALS_NUM_NODES_ERROR
                        , (TopologyData)
                      );

};

template <typename TopologyData, typename PermutationVector, unsigned NumNodes, unsigned Permutation>
struct check_permutations<TopologyData,PermutationVector,NumNodes,Permutation,Permutation>
{
  static const bool value = true;
};



template <typename TopologyData, unsigned Face = 0, unsigned NumFaces = TopologyData::num_faces >
struct check_faces
{

  typedef typename mpl::at_c<typename TopologyData::face_topology_vector, Face>::type face_topology_;
  typedef topology_data< face_topology_::value > face_topology;

  typedef typename mpl::at_c<typename TopologyData::face_node_ordinals_vector, Face>::type face_nodes;

  static const bool value =    (face_topology::num_nodes == mpl::size< face_nodes >::value)
                            && check_faces<TopologyData,Face+1>::value;

  BOOST_MPL_ASSERT_MSG(   value
                        , FACE_NODE_ORDINALS_NUM_NODES_ERROR
                        , (TopologyData)
                      );

};

template <typename TopologyData, unsigned Face>
struct check_faces<TopologyData,Face,Face>
{
  static const bool value = true;
};

template <typename TopologyData>
struct validate_topology_data
{
  static const bool edge_node_ordinals = TopologyData::num_edges == mpl::size<typename TopologyData::edge_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   edge_node_ordinals
                        , EDGE_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool face_topology = TopologyData::num_faces == mpl::size<typename TopologyData::face_topology_vector>::value;
  BOOST_MPL_ASSERT_MSG(   face_topology
                        , FACE_TOPOLOGY_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool face_node_ordinals = TopologyData::num_faces == mpl::size<typename TopologyData::face_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   face_node_ordinals
                        , FACE_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool permutation_node_ordinals = TopologyData::num_permutations == mpl::size<typename TopologyData::permutation_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   permutation_node_ordinals
                        , PERMUTATION_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool check_permutations_nodes = check_permutations<TopologyData, typename TopologyData::permutation_node_ordinals_vector, TopologyData::num_nodes>::value;

  static const bool check_face_nodes = check_faces<TopologyData>::value;

  static const bool value =    edge_node_ordinals
                            && face_topology
                            && face_node_ordinals
                            && permutation_node_ordinals
                            && check_permutations_nodes
                            && check_face_nodes
                           ;

};

} //unnamed namespace

STKUNIT_UNIT_TEST( stk_topology, validate_topology)
{
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::NODE          > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_2        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_3        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_3         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_4         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_6         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_4        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_8        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_9        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::PARTICLE      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_2_1D     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_3_1D     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::BEAM_2        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::BEAM_3        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_LINE_2  > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_LINE_3  > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_3_2D      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_4_2D      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_6_2D      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_4_2D     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_8_2D     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_9_2D     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_3   > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_4   > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_6   > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_4  > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_8  > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_9  > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_4         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_8         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_10        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_11        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_5     > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_13    > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_14    > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_6       > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_15      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_18      > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_8         > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_20        > >::value );
  STKUNIT_EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_27        > >::value );

}



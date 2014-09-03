#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <vector>

#ifndef __IBMCPP__

namespace {

TEST(side_node_ordinals, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;
  const unsigned nodes_per_side = 4;
  unsigned side1_node_ordinals[nodes_per_side] = {0, 0, 0, 0};

  unsigned expected_side1_node_ordinals[nodes_per_side] = {1, 2, 6, 5};

  //important usage note: side-ordinal is a 0-based number.
  //i.e., side 1 is labeled as side 2 in the exodus manual. it is the second side of the element.
  //note also that the side-node-ordinals returned are 0-based whereas the node labels in the
  //exodus manual diagrams are 1-based.
 
  const unsigned side1_ordinal = 1;

  //the following call to hex8.side_node_ordinals(..) fills our side1_node_ordinals array.
  hex8.side_node_ordinals(side1_ordinal, side1_node_ordinals);

  EXPECT_EQ(expected_side1_node_ordinals[0], side1_node_ordinals[0]);
  EXPECT_EQ(expected_side1_node_ordinals[1], side1_node_ordinals[1]);
  EXPECT_EQ(expected_side1_node_ordinals[2], side1_node_ordinals[2]);
  EXPECT_EQ(expected_side1_node_ordinals[3], side1_node_ordinals[3]);

  //stk::topology::HEX_8 (and other topologies) define tables of edge and face node ordinals using
  //boost::mpl::vector objects, which are compile-time meta-programming constructs. What this means
  //is that they are typedefs, not instance data. i.e., there is not an actual table of numbers sitting
  //in memory somewhere.
  //
  //The above hex8.side_node_ordinals method uses meta-programming (meta-functions, etc) to translate
  //the contents of the typedefs into values that are put into the user-supplied array.
  //
  //Below is a sample of what that nicely-hidden meta-programming looks like.
  //We will use boost::mpl 'random access iterators' to repeat the same EXPECT_EQ tests that
  //were performed above for the node-ordinals of side 1.
  //
  typedef typename stk::topology::topology_type<stk::topology::HEX_8>::face_node_ordinals_vector face_node_ordinals_vector;

  typedef boost::mpl::at_c<face_node_ordinals_vector, side1_ordinal>::type side1_node_ordinals_iterator;


  EXPECT_EQ( expected_side1_node_ordinals[0], (boost::mpl::at< side1_node_ordinals_iterator, boost::mpl::int_<0> >::type::value) );
  EXPECT_EQ( expected_side1_node_ordinals[1], (boost::mpl::at< side1_node_ordinals_iterator, boost::mpl::int_<1> >::type::value) );
  EXPECT_EQ( expected_side1_node_ordinals[2], (boost::mpl::at< side1_node_ordinals_iterator, boost::mpl::int_<2> >::type::value) );
  EXPECT_EQ( expected_side1_node_ordinals[3], (boost::mpl::at< side1_node_ordinals_iterator, boost::mpl::int_<3> >::type::value) );
}

}

#endif // __IBMCPP__

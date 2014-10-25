#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

TEST( exp_stk_topology, hex8_face_node_ordinals)
{
  unsigned expected_face_node_ordinals[] = { 1, 2, 6, 5};

  unsigned face_ordinal = 1;

  stk::topology hex8 = stk::topology::HEX_8;

  const unsigned nodes_per_face = 4;
  unsigned face_node_ordinals[nodes_per_face] = {0, 0, 0, 0};

  hex8.face_node_ordinals(face_ordinal, face_node_ordinals);

  EXPECT_EQ(expected_face_node_ordinals[0], face_node_ordinals[0]);
  EXPECT_EQ(expected_face_node_ordinals[1], face_node_ordinals[1]);
  EXPECT_EQ(expected_face_node_ordinals[2], face_node_ordinals[2]);
  EXPECT_EQ(expected_face_node_ordinals[3], face_node_ordinals[3]);
}


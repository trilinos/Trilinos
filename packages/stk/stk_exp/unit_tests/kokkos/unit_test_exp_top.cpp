#include <gtest/gtest.h>

struct Hex8 {
  static const unsigned faces_per_element = 6;
  static const unsigned nodes_per_face = 4;

  inline
  static const unsigned * face_node_ordinals(unsigned face_ordinal)
  {
    static const unsigned m_face_node_ordinals[faces_per_element][nodes_per_face] =
     { { 0, 1, 5, 4 },
       { 1, 2, 6, 5 },
       { 2, 3, 7, 6 },
       { 0, 4, 7, 3 },
       { 0, 3, 2, 1 },
       { 4, 5, 6, 7 } };
  
    return m_face_node_ordinals[face_ordinal];
  }
};

TEST( exp_topology, hex8_face_node_ordinals)
{
  unsigned expected_face_node_ordinals[] = { 1, 2, 6, 5};

  unsigned face_ordinal = 1;

  const unsigned* face_node_ordinals = Hex8::face_node_ordinals(face_ordinal);

  EXPECT_EQ(expected_face_node_ordinals[0], face_node_ordinals[0]);
  EXPECT_EQ(expected_face_node_ordinals[1], face_node_ordinals[1]);
  EXPECT_EQ(expected_face_node_ordinals[2], face_node_ordinals[2]);
  EXPECT_EQ(expected_face_node_ordinals[3], face_node_ordinals[3]);
}


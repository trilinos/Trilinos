// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QuadFace.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <array>
#include <unordered_set>
#include <vector>

namespace krino
{

//--------------------------------------------------------------------------
// QuadFace packs the four face-node local offsets into a single 128-bit key.
// These pin the value-type semantics that the hashing/comparison containers
// rely on -- no mesh required.
//--------------------------------------------------------------------------
TEST(QuadFace, defaultConstructed_isInvalid)
{
  QuadFace face;
  EXPECT_FALSE(face.is_valid());
  EXPECT_EQ(QuadFace::InvalidQuadFace, face.value());
}

TEST(QuadFace, nonzeroValue_isValid)
{
  QuadFace face(QuadFace::quad_face_value_type(7u));
  EXPECT_TRUE(face.is_valid());
}

TEST(QuadFace, equalityAndOrdering_followUnderlyingValue)
{
  QuadFace a(QuadFace::quad_face_value_type(5u));
  QuadFace b(QuadFace::quad_face_value_type(5u));
  QuadFace c(QuadFace::quad_face_value_type(7u));

  EXPECT_TRUE(a == b);
  EXPECT_TRUE(a != c);
  EXPECT_TRUE(a < c);
  EXPECT_FALSE(c < a);

  // The specialized std::hash forwards to the bitset hash, so equal faces hash
  // equally and can be de-duplicated in an unordered_set.
  std::unordered_set<QuadFace> faces;
  faces.insert(a);
  faces.insert(b); // same value as a
  faces.insert(c);
  EXPECT_EQ(2u, faces.size());
}

//--------------------------------------------------------------------------
// Mesh-based quad-face encoding/decoding, built on a single regular hex.
//--------------------------------------------------------------------------
class RegularHexQuadFace : public StkMeshHexFixture
{
protected:
  RegularHexQuadFace()
  {
    RegularHex meshSpec;
    set_valid_proc_sizes_for_test({1});
    if (is_valid_proc_size_for_test())
      StkMeshHexFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Entity get_element() { return get_owned_elements()[0]; }
};

// quad_face_from_unordered_nodes sorts its four inputs, so any permutation of
// the same four nodes yields the same face; get_quad_face_nodes_sorted_by_id
// decodes them back in id order.
TEST_F(RegularHexQuadFace, faceFromNodes_isPermutationInvariantAndRoundTrips)
{
  if (!is_valid_proc_size_for_test())
    return;

  const stk::mesh::BulkData & mesh = get_mesh();
  const stk::mesh::Entity * const nodes = mesh.begin_nodes(get_element());
  const std::array<stk::mesh::Entity, 4> faceNodes{nodes[0], nodes[1], nodes[2], nodes[3]};

  const QuadFace face = quad_face_from_unordered_nodes(mesh, faceNodes);
  EXPECT_TRUE(face.is_valid());

  // A different input ordering of the same four nodes must produce the same face.
  const QuadFace facePermuted = quad_face_from_unordered_nodes(mesh, nodes[2], nodes[0], nodes[3], nodes[1]);
  EXPECT_EQ(face, facePermuted);

  // Decoding recovers the four nodes sorted by identifier.
  const std::array<stk::mesh::Entity, 4> decoded = get_quad_face_nodes_sorted_by_id(face);
  for (unsigned i = 0; i + 1 < decoded.size(); ++i)
    EXPECT_LT(mesh.identifier(decoded[i]), mesh.identifier(decoded[i+1]));

  std::array<stk::mesh::EntityId, 4> expectedIds{ mesh.identifier(nodes[0]), mesh.identifier(nodes[1]),
                                                  mesh.identifier(nodes[2]), mesh.identifier(nodes[3]) };
  std::sort(expectedIds.begin(), expectedIds.end());
  for (unsigned i = 0; i < decoded.size(); ++i)
    EXPECT_EQ(expectedIds[i], mesh.identifier(decoded[i]));
}

// A hex has 6 quad faces; append/fill_entity_quad_faces must produce that many
// unique, valid faces.
TEST_F(RegularHexQuadFace, entityQuadFaces_areSixUniqueValidFacesForHex)
{
  if (!is_valid_proc_size_for_test())
    return;

  const stk::mesh::BulkData & mesh = get_mesh();

  std::vector<QuadFace> faces;
  fill_entity_quad_faces(mesh, get_element(), faces);
  ASSERT_EQ(6u, faces.size());
  for (const QuadFace & face : faces)
    EXPECT_TRUE(face.is_valid());

  const std::unordered_set<QuadFace> uniqueFaces(faces.begin(), faces.end());
  EXPECT_EQ(6u, uniqueFaces.size());
}

}

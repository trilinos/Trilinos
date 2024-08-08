// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <gtest/gtest.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/BuildMesh.hpp>

#include <unistd.h>

using stk::unit_test_util::build_mesh;

class StkTopologyTest : public stk::unit_test_util::MeshFixture
{
public:
  StkTopologyTest()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void init_mesh_with_wedge12_element(stk::mesh::BulkData& bulk)
  {
    std::string meshDesc = "0,1,WEDGE_12,1,2,3,4,5,6,7,8,9,10,11,12";
    stk::unit_test_util::setup_text_mesh(bulk, meshDesc);
    stk::mesh::create_all_sides(bulk, bulk.mesh_meta_data().universal_part(), {}, true);
  }

  stk::topology get_element_topology(const stk::mesh::BulkData& bulk)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(bulk, stk::topology::ELEMENT_RANK, elements);
    EXPECT_EQ(elements.size(), 1u);
    return bulk.bucket(elements[0]).topology();
  }

  void output_mesh(stk::mesh::BulkData& bulk, const std::string& fileName)
  {
    stk::io::write_mesh(fileName, bulk);
  }

  void check_mesh_has_wedge12(stk::mesh::BulkData& bulk)
  {
    EXPECT_EQ(get_element_topology(bulk), stk::topology::WEDGE_12);

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    EXPECT_TRUE(bulk.is_valid(elem));
    const stk::mesh::Entity* faces = bulk.begin_faces(elem);
    EXPECT_EQ(5u, bulk.num_faces(elem));
    EXPECT_EQ(bulk.bucket(faces[0]).topology(), stk::topology::QUAD_6);
    EXPECT_EQ(bulk.bucket(faces[1]).topology(), stk::topology::QUAD_6);
    EXPECT_EQ(bulk.bucket(faces[2]).topology(), stk::topology::QUAD_6);
    EXPECT_EQ(bulk.bucket(faces[3]).topology(), stk::topology::TRI_6);
    EXPECT_EQ(bulk.bucket(faces[4]).topology(), stk::topology::TRI_6);
  }
};

TEST_F(StkTopologyTest, createWedge12OnMesh)
{
  if(get_parallel_size() != 1) { return; }

  init_mesh_with_wedge12_element(get_bulk());
  check_mesh_has_wedge12(get_bulk());
}

TEST_F(StkTopologyTest, writeWedge12ToExodusFile)
{
  if(get_parallel_size() != 1) { return; }

  const std::string fileName = "single_wedge_case.e";
  init_mesh_with_wedge12_element(get_bulk());

  output_mesh(get_bulk(), fileName);

  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(3, get_comm());
  stk::mesh::BulkData& newBulk = *newBulkPtr;
  stk::mesh::MetaData& newMeta = newBulk.mesh_meta_data();
  stk::io::fill_mesh(fileName, newBulk);
  stk::mesh::create_all_sides(newBulk, newMeta.universal_part(), {}, true);

  check_mesh_has_wedge12(newBulk);

  unlink(fileName.c_str());
}

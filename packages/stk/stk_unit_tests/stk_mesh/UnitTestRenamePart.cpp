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
// 

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Part.hpp>       // for PartLess, intersect, Part, etc
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <string>                       // for string, char_traits

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::impl::PartRepository;

namespace {

TEST(RenamePart, same_name)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  Part & myPart = m.declare_part_with_topology("myPart", stk::topology::HEX_8);
  m.commit();

  EXPECT_NO_THROW(m.rename(myPart, myPart.name()));
}

TEST(RenamePart, internal_part_throws)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  m.commit();

  EXPECT_ANY_THROW(m.rename(m.universal_part(), "foo"));
}

TEST(RenamePart, dont_use_internal_part_name)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  Part & myPart = m.declare_part_with_topology("myPart", stk::topology::HEX_8);
  m.commit();

  std::string newName("{illegalName}");
  EXPECT_ANY_THROW(m.rename(myPart, newName));
}

TEST(RenamePart, dont_use_name_of_other_part)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  Part & myPart = m.declare_part_with_topology("myPart", stk::topology::HEX_8);
  Part & myOtherPart = m.declare_part_with_topology("myOtherPart", stk::topology::HEX_8);
  m.commit();

  const std::string& newName = myOtherPart.name();
  EXPECT_ANY_THROW(m.rename(myPart, newName));
}

TEST(RenamePart, not_during_mesh_mod)
{
  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(stk::parallel_machine_world())
                                                      .set_spatial_dimension(spatial_dimension)
                                                      .create();
  stk::mesh::MetaData& m = bulkPtr->mesh_meta_data();
  Part & myPart = m.declare_part_with_topology("myPart", stk::topology::HEX_8);
  m.commit();

  bulkPtr->modification_begin();
  const std::string newName("foo");
  EXPECT_ANY_THROW(m.rename(myPart, newName));
}

#ifndef NDEBUG
TEST(RenamePart, name_not_parallel_consistent_throws)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 2) { GTEST_SKIP(); }

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(stk::parallel_machine_world())
                                                      .set_spatial_dimension(spatial_dimension)
                                                      .create();
  stk::mesh::MetaData& m = bulkPtr->mesh_meta_data();
  Part & myPart = m.declare_part_with_topology("myPart", stk::topology::HEX_8);
  m.commit();

  const int myProc = stk::parallel_machine_rank(stk::parallel_machine_world());
  const std::string newName = myProc==0 ? "newPartName" : "differentName";
  EXPECT_ANY_THROW(m.rename(myPart, newName));
}
#endif

TEST(RenamePart, get_same_part_with_new_name)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  std::string origName("myPart");
  Part & myPart = m.declare_part_with_topology(origName, stk::topology::HEX_8);
  m.commit();

  const std::string newName("newPartName");
  EXPECT_NO_THROW(m.rename(myPart, newName));

  Part* partWithNewName = m.get_part(newName);
  ASSERT_TRUE(partWithNewName != nullptr);
  EXPECT_EQ(partWithNewName->mesh_meta_data_ordinal(), myPart.mesh_meta_data_ordinal());
  EXPECT_EQ(newName, partWithNewName->name());

  Part* partWithOrigName = m.get_part(origName);
  EXPECT_EQ(nullptr, partWithOrigName);
}

TEST(RenamePart, rename_part_with_alias)
{
  const int spatial_dimension = 3;
  stk::mesh::MetaData m(spatial_dimension);
  std::string origName("myPart");
  Part & myPart = m.declare_part_with_topology(origName, stk::topology::HEX_8);
  m.add_part_alias(myPart, origName);
  m.commit();

  const std::string newName("newPartName");
  EXPECT_NO_THROW(m.rename(myPart, newName));

  Part* partWithOrigName = m.get_part(origName);
  EXPECT_EQ(nullptr, partWithOrigName);
}

} // anonymous namespace


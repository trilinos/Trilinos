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

#include <gtest/gtest.h>                // for AssertHelper, TEST, etc
#include <stddef.h>                     // for NULL, size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error, logic_error
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <string>                       // for string, operator==
#include <vector>                       // for vector
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/FillMesh.hpp"

namespace stk { namespace mesh { class Part; } }

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::MeshBuilder;

namespace {

TEST( UnitTestPartAlias, noAliasForDefaultCreatedPart )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  EXPECT_TRUE(aliases.empty());
}

TEST( UnitTestPartAlias, getAliasForAliasedPart )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias = "my_part";
  metadata.add_part_alias(pa, alias);

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(1u, aliases.size());
  EXPECT_EQ(alias, aliases[0]);
}

TEST( UnitTestPartAlias, partAliasesAreUnique )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );
  Part &pb = metadata.declare_part( std::string("b") , node_rank );

  const std::string alias = "my_part";
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias));
  EXPECT_THROW(metadata.add_part_alias(pb, alias), std::logic_error);
}

TEST( UnitTestPartAlias, partAliasReregistration )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias = "my_part";
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias));
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias));

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(1u, aliases.size());
  EXPECT_EQ(alias, aliases[0]);
}

TEST( UnitTestPartAlias, multiplePartAliasRegistration )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias1 = "alias_1";
  const std::string alias2 = "alias_2";
  const std::string alias3 = "alias_3";
  const std::string alias4 = "alias_4";
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias1));
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias2));
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias3));
  EXPECT_NO_THROW(metadata.add_part_alias(pa, alias4));

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(4u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);
  EXPECT_EQ(alias2, aliases[1]);
  EXPECT_EQ(alias3, aliases[2]);
  EXPECT_EQ(alias4, aliases[3]);
}

TEST( UnitTestPartAlias, partAliasesAreCaseSensitive )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias1 = "MY_PART";
  const std::string alias2 = "my_part";
  metadata.add_part_alias(pa, alias1);
  metadata.add_part_alias(pa, alias2);

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(2u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);
  EXPECT_EQ(alias2, aliases[1]);
}

TEST( UnitTestPartAlias, deletePartAliasExact )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias1 = "MY_PART";
  const std::string alias2 = "MY_part";
  const std::string alias3 = "my_part";
  metadata.add_part_alias(pa, alias1);
  metadata.add_part_alias(pa, alias2);
  metadata.add_part_alias(pa, alias3);

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(3u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);
  EXPECT_EQ(alias2, aliases[1]);
  EXPECT_EQ(alias3, aliases[2]);

  EXPECT_NO_THROW(metadata.delete_part_alias(pa, alias2));

  aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(2u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);
  EXPECT_EQ(alias3, aliases[1]);

  EXPECT_NO_THROW(metadata.delete_part_alias(pa, alias3));

  aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(1u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);

  EXPECT_NO_THROW(metadata.delete_part_alias(pa, alias1));

  aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(0u, aliases.size());
}

TEST( UnitTestPartAlias, deletePartAliasCaseInsensitive )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );

  const std::string alias1 = "MY_PART";
  const std::string alias2 = "MY_part";
  const std::string alias3 = "my_part";
  metadata.add_part_alias(pa, alias1);
  metadata.add_part_alias(pa, alias2);
  metadata.add_part_alias(pa, alias3);

  std::vector<std::string> aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(3u, aliases.size());
  EXPECT_EQ(alias1, aliases[0]);
  EXPECT_EQ(alias2, aliases[1]);
  EXPECT_EQ(alias3, aliases[2]);

  EXPECT_NO_THROW(metadata.delete_part_alias_case_insensitive(pa, alias2));

  aliases = metadata.get_part_aliases(pa);
  ASSERT_EQ(0u, aliases.size());
}

TEST( UnitTestPartAlias, getPartFromAlias )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  std::string partName_a("a");
  std::string partName_b("b");
  Part &pa = metadata.declare_part( partName_a , node_rank );
  Part &pb = metadata.declare_part( partName_b , node_rank );

  const std::string alias1 = "MY_PART";
  const std::string alias2 = "MY_part";
  metadata.add_part_alias(pa, alias1);
  metadata.add_part_alias(pa, alias2);

  const std::string alias3 = "my_part";
  metadata.add_part_alias(pb, alias3);

  Part* part_a       = metadata.get_part(partName_a);
  Part* aliasedPart1 = metadata.get_part(alias1);
  Part* aliasedPart2 = metadata.get_part(alias2);

  Part* part_b       = metadata.get_part(partName_b);
  Part* aliasedPart3 = metadata.get_part(alias3);

  ASSERT_TRUE(nullptr != part_a);
  ASSERT_TRUE(nullptr != aliasedPart1);
  ASSERT_TRUE(nullptr != aliasedPart2);

  ASSERT_TRUE(nullptr != part_b);
  ASSERT_TRUE(nullptr != aliasedPart3);

  EXPECT_EQ(pa.mesh_meta_data_ordinal(), part_a->mesh_meta_data_ordinal());
  EXPECT_EQ(pa.mesh_meta_data_ordinal(), aliasedPart1->mesh_meta_data_ordinal());
  EXPECT_EQ(pa.mesh_meta_data_ordinal(), aliasedPart2->mesh_meta_data_ordinal());

  EXPECT_EQ(pb.mesh_meta_data_ordinal(), part_b->mesh_meta_data_ordinal());
  EXPECT_EQ(pb.mesh_meta_data_ordinal(), aliasedPart3->mesh_meta_data_ordinal());
}

}
//----------------------------------------------------------------------


/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for NULL, size_t
#include <exception>                    // for exception
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <string>                       // for string, operator==, etc
#include <vector>                       // for vector
#include "Shards_CellTopologyData.h"    // for CellTopologyData
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/NamedPair.hpp"
namespace stk { namespace mesh { class Part; } }






using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using std::cout;
using std::endl;

//----------------------------------------------------------------------

namespace {

TEST( UnitTestMetaData, testMetaData )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  MetaData metadata_committed(spatial_dimension);
  MetaData metadata_not_committed(spatial_dimension);
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );
  Part &pb = metadata.declare_part( std::string("b") , node_rank );
  Part &pc = metadata.declare_part( std::string("c") , node_rank );
  Part &pd = metadata.declare_part( std::string("d") , node_rank );
  Part &pe = metadata.declare_part( std::string("e") , node_rank );
  PartVector part_vector;
  metadata_committed.commit();

  //test get_part with part that does not exist
  std::string test_string = "this_part_does_not_exist";
  ASSERT_THROW( metadata_committed.get_part(test_string,"test_throw"),std::runtime_error);

  //test get_part with valid part
  ASSERT_TRUE( metadata.get_part(std::string("a"),"do_not_throw"));



  part_vector.push_back(& pa);
  part_vector.push_back(& pb);
  part_vector.push_back(& pc);
  part_vector.push_back(& pd);

  //Test declare_part_subset
  ASSERT_THROW(  metadata.declare_part_subset( pe, pe), std::runtime_error);

  metadata.commit();
}

TEST( UnitTestMetaData, rankHigherThanDefined )
{
  //Test function entity_rank_name in MetaData.cpp
  const int spatial_dimension = 3;
  const std::vector<std::string> & rank_names = stk::mesh::entity_rank_names();
  MetaData metadata(spatial_dimension, rank_names);

  const std::string& i_name2 =  metadata.entity_rank_name( stk::topology::EDGE_RANK );

  ASSERT_TRUE( i_name2 == rank_names[stk::topology::EDGE_RANK] );

  EntityRank one_rank_higher_than_defined = static_cast<EntityRank>(rank_names.size());

  ASSERT_THROW(
    metadata.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
                        );
}

TEST( UnitTestMetaData, testEntityRepository )
{
  static const size_t spatial_dimension = 3;

  //Test Entity repository - covering EntityRepository.cpp/hpp
  stk::mesh::MetaData meta ( spatial_dimension );
  stk::mesh::Part & part = meta.declare_part("another part");
  stk::mesh::Part & hex_part = meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);

  meta.commit();

  stk::mesh::BulkData bulk ( meta , MPI_COMM_WORLD );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &part );
  std::vector<stk::mesh::Part *> elem_parts;
  elem_parts.push_back( &part );
  elem_parts.push_back( &hex_part );

  int rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  int size = stk::parallel_machine_size( MPI_COMM_WORLD );
  PartVector tmp(1);

  bulk.modification_begin();

  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::Entity node = stk::mesh::Entity();
  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    node = bulk.declare_entity( stk::topology::NODE_RANK , new_id+1 , add_part );
    nodes.push_back(node);
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity elem  = bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id+1 , elem_parts );

  for (unsigned ord = 0; ord < 8; ++ord)
  {
    bulk.declare_relation(elem, nodes[ord], ord);
  }

  bulk.entity_comm_map_clear(bulk.entity_key(elem));

  bulk.entity_comm_map_clear_ghosting(bulk.entity_key(elem));

  const stk::mesh::Ghosting & ghost = bulk.aura_ghosting();

  bulk.modification_end();

  ASSERT_FALSE(bulk.entity_comm_map_erase(bulk.entity_key(elem), ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  ASSERT_FALSE(bulk.entity_comm_map_erase(bulk.entity_key(elem), comm_info));
  ASSERT_TRUE(bulk.entity_comm_map_insert(elem, comm_info));
}

TEST( UnitTestMetaData, noEntityTypes )
{
  //MetaData constructor fails because there are no entity types:
  std::vector<std::string> wrong_names(1, "foo");
  ASSERT_THROW(
    MetaData metadata(3 /*dim*/, wrong_names),
    std::runtime_error
    );
}
TEST( UnitTestMetaData, declare_part_with_rank )
{
  //MetaData constructor fails because there are no entity types:
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  metadata.declare_part("foo");
  ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));
  ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));

  // Should throw because we're trying to change rank
  ASSERT_THROW(metadata.declare_part("foo",stk::topology::FACE_RANK),std::runtime_error);

  // Should not throw since we did not provide rank
  metadata.declare_part("foo");
}

TEST( UnitTestMetaData, declare_attribute_no_delete )
{
  //Coverage of declare_attribute_no_delete in MetaData.hpp
  const CellTopologyData * singleton = NULL;
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  Part &pa = metadata.declare_part( std::string("a") , stk::topology::NODE_RANK );
  metadata.declare_attribute_no_delete( pa, singleton);
  metadata.commit();
}

}
//----------------------------------------------------------------------






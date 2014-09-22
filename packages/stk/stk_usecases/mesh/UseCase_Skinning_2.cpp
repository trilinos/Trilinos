// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <iostream>

#include <mesh/UseCase_Skinning.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

unsigned count_skin_entities( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, stk::mesh::EntityRank skin_rank)
{
  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(mesh);

  stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const stk::mesh::BucketVector& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

// \TODO This function is a general utility, to be moved into stk_mesh and unit tested.

// Destroy this entity and any lower ranking entities in this entity's closure that no longer
// have any upward relations.
void destroy_entity_closure( stk::mesh::BulkData & mesh, stk::mesh::Entity entity)
{
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);

  for (stk::mesh::EntityRank irank = stk::topology::END_RANK;
         irank != stk::topology::BEGIN_RANK;)
  {
    --irank;

    stk::mesh::Entity const *relations = mesh.begin(entity, irank);
    int num_relations = mesh.num_connectivity(entity, irank);
    stk::mesh::ConnectivityOrdinal const *relation_ordinals = mesh.begin_ordinals(entity, irank);

    ThrowErrorMsgIf( (num_relations > 0) && (irank > entity_rank),
                     "Unable to destroy and entity with upward relations" );

    for (int j = num_relations - 1; j >= 0; --j)
    {
      stk::mesh::Entity related_entity = relations[j];

      if (!related_entity.is_local_offset_valid())
        continue;

      stk::mesh::RelationIdentifier rel_id = relation_ordinals[j];
      stk::mesh::EntityRank related_entity_rank = irank;

      mesh.destroy_relation( entity, related_entity, rel_id );

      bool related_entity_no_upward_rels = true;
      for (stk::mesh::EntityRank krank = ++related_entity_rank;
            krank != stk::topology::END_RANK;
            ++krank)
      {
        if (mesh.count_valid_connectivity(related_entity, krank) > 0)
        {
          related_entity_no_upward_rels = false;
          break;
        }
      }

      if (related_entity_no_upward_rels)
      {
        destroy_entity_closure(mesh, related_entity);
      }
    }
  }

  mesh.destroy_entity(entity);
}

}

// \TODO ASCII art to illustrate the whole use case geometry
//       both before and after.

bool skinning_use_case_2(stk::ParallelMachine pm)
{
  const unsigned nx = 2 , ny = 1 , nz = 1 ;

  bool result = true;

  //TODO check the skin after each update to ensure that the appropriate
  //number of faces  and particles exist.
  try {
    stk::mesh::fixtures::HexFixture fixture( pm , nx , ny , nz );
    const stk::mesh::EntityRank side_rank = fixture.m_meta.side_rank();

    const int p_rank = fixture.m_bulk_data.parallel_rank();
    const int p_size = fixture.m_bulk_data.parallel_size();

    stk::mesh::Part & skin_part = fixture.m_meta.declare_part("skin_part");

    stk::mesh::Part & shell_part = fixture.m_meta.declare_part_with_topology("shell_part", stk::topology::SHELL_QUAD_4);

    fixture.m_meta.commit();

    fixture.generate_mesh();

    fixture.m_bulk_data.modification_begin();

    if ( p_rank + 1 == p_size ) {
      // Verifies when all three elements on different processes, for p_size > 2
      //add shell between the two elements

      stk::mesh::EntityId elem_node[4] ;

      // Query nodes from this simple grid fixture via the (i,j,k) indices.
      elem_node[0] = fixture.node_id( 1, 0, 0 );
      elem_node[1] = fixture.node_id( 1, 1, 0 );
      elem_node[2] = fixture.node_id( 1, 1, 1 );
      elem_node[3] = fixture.node_id( 1, 0, 1 );

      stk::mesh::EntityId elem_id = 3;

      stk::mesh::declare_element( fixture.m_bulk_data, shell_part, elem_id, elem_node);
    }
    fixture.m_bulk_data.modification_end();

    {
      stk::mesh::PartVector add_parts(1,&skin_part);
      stk::mesh::skin_mesh(fixture.m_bulk_data, add_parts);
    }

    //----------------------------------------------------------------------
    //Actual usecase
    //----------------------------------------------------------------------

    {
      int num_skin_entities = count_skin_entities( fixture.m_bulk_data, skin_part, side_rank);

      stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

      if ( num_skin_entities != 10 ) {
        result = false;
        std::cerr << std::endl << "incorrect number of entities in skin.  Expected 10, Found "
          << num_skin_entities << std::endl;
      }
    }

    // Kill element on the "left" of the shell:
    fixture.m_bulk_data.modification_begin();
    stk::mesh::Entity elem_to_kill = fixture.elem( 0 , 0 , 0 ); // (i,j,k) indices
    if ( fixture.m_bulk_data.is_valid(elem_to_kill) && p_rank == fixture.m_bulk_data.parallel_owner_rank(elem_to_kill) ) {
      // Destroy element and its sides and nodes
      // that are not in the closure of another element.
      destroy_entity_closure( fixture.m_bulk_data, elem_to_kill);
    }

    fixture.m_bulk_data.modification_end();

    {
      stk::mesh::PartVector add_parts(1,&skin_part);
      stk::mesh::skin_mesh(fixture.m_bulk_data, add_parts);
    }

    {
      int num_skin_entities = count_skin_entities( fixture.m_bulk_data, skin_part, side_rank);

      stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

      // Verify that the correct 6 sides are present.

      if ( num_skin_entities != 6 ) {
        result = false;
        std::cerr << std::endl << "incorrect number of entities in skin.  Expected 6, Found "
          << num_skin_entities << std::endl;
      }
    }
  }
  catch(std::exception & e) {
    std::cerr << std::endl << e.what() << std::endl;
    result = false;
  }

  return result;
}

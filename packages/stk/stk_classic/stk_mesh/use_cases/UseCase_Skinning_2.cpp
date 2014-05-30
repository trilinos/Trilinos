#include <iostream>

#include <use_cases/UseCase_Skinning.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

unsigned count_skin_entities( stk_classic::mesh::BulkData & mesh, stk_classic::mesh::Part & skin_part, stk_classic::mesh::EntityRank skin_rank)
{
  const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(mesh);

  stk_classic::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk_classic::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

// \TODO This function is a general utility, to be moved into stk_mesh and unit tested.

// Destroy this entity and any lower ranking entities in this entity's closure that no longer
// have any upward relations.
void destroy_entity_closure( stk_classic::mesh::BulkData & mesh, stk_classic::mesh::Entity * entity)
{
  stk_classic::mesh::PairIterRelation relations = entity->relations();
  stk_classic::mesh::EntityRank entity_rank = entity->entity_rank();

  ThrowErrorMsgIf( !relations.empty() &&
                   relations.back().entity()->entity_rank() > entity_rank,
                   "Unable to destroy and entity with upward relations" );

  for (; !entity->relations().empty();) {
    stk_classic::mesh::Entity * related_entity = (entity->relations().back().entity());
    stk_classic::mesh::RelationIdentifier rel_id = entity->relations().back().identifier();
    stk_classic::mesh::EntityRank related_entity_rank = related_entity->entity_rank();

    mesh.destroy_relation( *entity, *related_entity, rel_id );

    stk_classic::mesh::PairIterRelation related_entity_relations = related_entity->relations();

    //  Only destroy if there are no upward relations
    if ( related_entity_relations.empty() ||
        related_entity_relations.back().entity()->entity_rank() < related_entity_rank )
    {
      destroy_entity_closure(mesh,related_entity);
    }
  }

  mesh.destroy_entity(entity);
}

}

// \TODO ASCII art to illustrate the whole use case geometry
//       both before and after.

bool skinning_use_case_2(stk_classic::ParallelMachine pm)
{
  const unsigned nx = 2 , ny = 1 , nz = 1 ;

  bool result = true;

  //TODO check the skin after each update to ensure that the appropriate
  //number of faces  and particles exist.
  try {
    stk_classic::mesh::fixtures::HexFixture fixture( pm , nx , ny , nz );
    const stk_classic::mesh::EntityRank element_rank = fixture.m_fem_meta.element_rank();
    const stk_classic::mesh::EntityRank side_rank = fixture.m_fem_meta.side_rank();

    const unsigned p_rank = fixture.m_bulk_data.parallel_rank();
    const unsigned p_size = fixture.m_bulk_data.parallel_size();

    stk_classic::mesh::Part & skin_part = fixture.m_fem_meta.declare_part("skin_part");

    stk_classic::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
    stk_classic::mesh::Part & shell_part = fixture.m_fem_meta.declare_part("shell_part", shell_top);

    fixture.m_fem_meta.commit();

    fixture.generate_mesh();

    fixture.m_bulk_data.modification_begin();

    if ( p_rank + 1 == p_size ) {
      // Verifies when all three elements on different processes, for p_size > 2
      //add shell between the two elements

      stk_classic::mesh::EntityId elem_node[4] ;

      // Query nodes from this simple grid fixture via the (i,j,k) indices.
      elem_node[0] = fixture.node_id( 1, 0, 0 );
      elem_node[1] = fixture.node_id( 1, 1, 0 );
      elem_node[2] = fixture.node_id( 1, 1, 1 );
      elem_node[3] = fixture.node_id( 1, 0, 1 );

      stk_classic::mesh::EntityId elem_id = 3;

      stk_classic::mesh::fem::declare_element( fixture.m_bulk_data, shell_part, elem_id, elem_node);
    }
    fixture.m_bulk_data.modification_end();

    stk_classic::mesh::skin_mesh(fixture.m_bulk_data, element_rank, &skin_part);

    //----------------------------------------------------------------------
    //Actual usecase
    //----------------------------------------------------------------------

    {
      int num_skin_entities = count_skin_entities( fixture.m_bulk_data, skin_part, side_rank);

      stk_classic::all_reduce(pm, stk_classic::ReduceSum<1>(&num_skin_entities));

      if ( num_skin_entities != 10 ) {
        result = false;
        std::cerr << std::endl << "incorrect number of entities in skin.  Expected 10, Found "
          << num_skin_entities << std::endl;
      }
    }

    // Kill element on the "left" of the shell:
    fixture.m_bulk_data.modification_begin();
    stk_classic::mesh::Entity * elem_to_kill = fixture.elem( 0 , 0 , 0 ); // (i,j,k) indices
    if ( elem_to_kill != NULL && p_rank == elem_to_kill->owner_rank() ) {
      // Destroy element and its sides and nodes
      // that are not in the closure of another element.
      destroy_entity_closure( fixture.m_bulk_data, elem_to_kill);
    }

    fixture.m_bulk_data.modification_end();

    stk_classic::mesh::skin_mesh( fixture.m_bulk_data, element_rank, &skin_part);

    {
      int num_skin_entities = count_skin_entities( fixture.m_bulk_data, skin_part, side_rank);

      stk_classic::all_reduce(pm, stk_classic::ReduceSum<1>(&num_skin_entities));

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

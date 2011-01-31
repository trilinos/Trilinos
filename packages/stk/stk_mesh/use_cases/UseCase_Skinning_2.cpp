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

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

unsigned count_skin_entities( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, stk::mesh::EntityRank skin_rank)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

// \TODO This function is a general utility, to be moved into stk_mesh and unit tested.

// Destroy this entity and any lower ranking entities in this entity's closure that no longer
// have any upward relations.
void destroy_entity_closure( stk::mesh::BulkData & mesh, stk::mesh::Entity * entity)
{
  stk::mesh::PairIterRelation relations = entity->relations();
  stk::mesh::EntityRank entity_rank = entity->entity_rank();

  ThrowErrorMsgIf( !relations.empty() &&
                   relations.back().entity()->entity_rank() > entity_rank,
                   "Unable to destroy and entity with upward relations" );

  for (; !entity->relations().empty();) {
    stk::mesh::Entity * related_entity = (entity->relations().back().entity());
    stk::mesh::EntityRank related_entity_rank = related_entity->entity_rank();

    mesh.destroy_relation( *entity, *related_entity);

    stk::mesh::PairIterRelation related_entity_relations = related_entity->relations();

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

bool skinning_use_case_2(stk::ParallelMachine pm)
{
  const unsigned nx = 2 , ny = 1 , nz = 1 ;

  bool result = true;

  //TODO check the skin after each update to ensure that the appropriate
  //number of faces  and particles exist.
  try {
    stk::mesh::fixtures::HexFixture fixture( pm , nx , ny , nz );
    const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fixture.m_fem);
    const stk::mesh::EntityRank side_rank = stk::mesh::fem::side_rank(fixture.m_fem);

    const unsigned p_rank = fixture.m_bulk_data.parallel_rank();
    const unsigned p_size = fixture.m_bulk_data.parallel_size();

    stk::mesh::Part & skin_part = declare_part(fixture.m_meta_data, "skin_part");

    stk::mesh::Part & shell_part = stk::mesh::declare_part<shards::ShellQuadrilateral<4> >(fixture.m_meta_data, "shell_part");

    fixture.m_meta_data.commit();

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

    stk::mesh::skin_mesh(fixture.m_bulk_data, element_rank, &skin_part);

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
    stk::mesh::Entity * elem_to_kill = fixture.elem( 0 , 0 , 0 ); // (i,j,k) indices
    if ( elem_to_kill != NULL && p_rank == elem_to_kill->owner_rank() ) {
      // Destroy element and its sides and nodes
      // that are not in the closure of another element.
      destroy_entity_closure( fixture.m_bulk_data, elem_to_kill);
    }

    fixture.m_bulk_data.modification_end();

    stk::mesh::skin_mesh( fixture.m_bulk_data, element_rank, &skin_part);

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

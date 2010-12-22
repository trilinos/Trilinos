/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Skinning.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stdexcept>
#include <iostream>

namespace {

void destroy_entity_and_create_particles(
    stk::mesh::fixtures::HexFixture & fixture,
    stk::mesh::Part & skin_part,
    stk::mesh::Entity * elem
    )
{
  stk::mesh::fem::FEMInterface &fem = fixture.m_fem;
  const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fem);

  const unsigned p_rank = fixture.m_bulk_data.parallel_rank();

  fixture.m_bulk_data.modification_begin();
  const stk::mesh::EntityRank particle_rank = element_rank;

  // forumlate request for 8 particles on owning process
  std::vector<size_t> requests(fixture.m_meta_data.entity_rank_count(), 0);
  if ( elem != NULL && p_rank == elem->owner_rank() ) {
    requests[particle_rank] = 8;
  }

  // create the particles
  stk::mesh::EntityVector new_particles;
  fixture.m_bulk_data.generate_new_entities(requests, new_particles);

  if ( ! new_particles.empty() ) {
    // Get node relations
    stk::mesh::PairIterRelation relations = elem->relations();

    std::vector<stk::mesh::Part*> add_parts;
    add_parts.push_back(&skin_part);

    // iterate over the new particles
    for (unsigned i = 0; i != 8; ++i) {
      // add the particles to the skin_part
      fixture.m_bulk_data.change_entity_parts( *new_particles[i],
                                             add_parts );
      // copy fields from nodes to particles
      fixture.m_bulk_data.copy_entity_fields( *(relations[i].entity()),
                                            *new_particles[i] );
    }
  }

  // delete element and entities in closure that have been orphaned
  if ( elem != NULL ) {
    stk::mesh::PairIterRelation relations = elem->relations();
    stk::mesh::EntityVector downward_relations;

    for (; relations.first != relations.second;) {
      --relations.second;
      stk::mesh::Entity * current_entity = (relations.second->entity());

      downward_relations.push_back(current_entity);
    }

    fixture.m_bulk_data.destroy_entity( elem );

    // destroy the related entities if they are not connected to any
    // higher order entities
    for (stk::mesh::EntityVector::iterator itr = downward_relations.begin();
        itr != downward_relations.end(); ++itr) {
      stk::mesh::Entity * current_entity = *itr;

      if (current_entity->relations(element_rank).empty()) {
        fixture.m_bulk_data.destroy_entity( current_entity );
      }
    }
  }

  fixture.m_bulk_data.modification_end();

  stk::mesh::skin_mesh( fixture.m_bulk_data, element_rank, &skin_part);
}

}

bool skinning_use_case_1b(stk::ParallelMachine pm)
{
  const unsigned nx = 3 , ny = 3 , nz = 3 ;

  bool result = true;

  //TODO check the skin after each update to ensure that the appropriate
  //number of faces  and particles exist.
  try {
    for ( unsigned iz = 0 ; iz < nz ; ++iz ) {
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
      stk::mesh::fixtures::HexFixture fixture( pm , nx , ny , nz );
      const stk::mesh::EntityRank element_rank = stk::mesh::fem::element_rank(fixture.m_fem);

      const stk::mesh::EntityRank particle_rank = element_rank;

      stk::mesh::Part & skin_part =
        declare_part(fixture.m_meta_data, "skin_part");

      stk::mesh::put_field( fixture.m_coord_field,
                            particle_rank,
                            fixture.m_meta_data.universal_part(),
                            3 );

      fixture.m_meta_data.commit();

      fixture.generate_mesh();

      stk::mesh::skin_mesh(fixture.m_bulk_data, element_rank, &skin_part);

      stk::mesh::Entity * elem = fixture.elem(ix , iy , iz);

      destroy_entity_and_create_particles(fixture, skin_part, elem);
    }
    }
    }
  }
  catch(std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
    result = false;
  }

  return result;
}

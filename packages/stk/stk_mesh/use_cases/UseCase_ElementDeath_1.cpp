/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/GridFixture.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

bool element_death_use_case(stk::ParallelMachine pm)
{

  GridFixture fixture(pm);
  stk::mesh::BulkData& bulk_data = fixture.bulk_data();
  stk::mesh::MetaData& meta_data = fixture.meta_data();

  // Arbitrarily chose faces 9, 10 to kill; see GridFixture.cpp
  std::vector<unsigned> entity_ids_to_kill;
  entity_ids_to_kill.push_back(9);
  entity_ids_to_kill.push_back(10);

  std::vector<stk::mesh::Entity*> entities_to_kill;
  for (std::vector<unsigned>::const_iterator itr = entity_ids_to_kill.begin();
         itr != entity_ids_to_kill.end(); ++itr) {
    stk::mesh::Entity * temp = bulk_data.get_entity(stk::mesh::Face, *itr);
    if (temp != NULL && temp->owner_rank() == bulk_data.parallel_rank()) {
      entities_to_kill.push_back(temp);
    }
  }

  // find the parallel-consistent closure of the elements to be killed
  std::vector<stk::mesh::Entity*> entities_closure;
  stk::mesh::find_closure(bulk_data,
                          entities_to_kill,
                          entities_closure);

  // find the boundary of the elements we're killing
  stk::mesh::EntitySideVector boundary;
  stk::mesh::boundary_analysis(bulk_data, entities_closure,
                               stk::mesh::Face, boundary);

  // Now do the element death. Kill entities by moving them to the dead part.
  bulk_data.modification_begin();
  std::vector<stk::mesh::Part*> add_parts;
  add_parts.push_back(fixture.dead_part());
  for (std::vector<stk::mesh::Entity*>::iterator itr = entities_to_kill.begin();
       itr != entities_to_kill.end(); ++itr) {
    bulk_data.change_entity_parts(**itr, add_parts);
  }

  // Counting the boundary sides
  unsigned new_side_count = 0;
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
       itr != boundary.end(); ++itr) {
    const stk::mesh::Entity* live_element = itr->outside.entity;
    if (live_element &&
        live_element->owner_rank() == bulk_data.parallel_rank()) {
      ++new_side_count;
    }
  }

  // Ask for new globally unique side ids
  std::vector<size_t> requests(meta_data.entity_type_count(), 0);
  std::vector<stk::mesh::Entity *> requested_entities;
  requests[stk::mesh::Edge] = new_side_count;
  bulk_data.generate_new_entities(requests, requested_entities);

  // Create boundaries between live and dead entities
  new_side_count = 0;
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
       itr != boundary.end(); ++itr) {
    stk::mesh::Entity* live_element = itr->outside.entity;
    const unsigned local_id = itr->outside.side_id;
    if (live_element &&
        live_element->owner_rank() == bulk_data.parallel_rank()) {
      stk::mesh::Entity & side = *requested_entities[new_side_count];
      stk::mesh::declare_element_side(*live_element, side,
                                      local_id,
                                      fixture.boundary_part());
      ++new_side_count;
    }
  }

  bulk_data.modification_end();

  // Check to see that results are correct
  // print live, print dead etc
  int passed = 1;
  stk::mesh::Selector selector = !(*fixture.dead_part()) & (*fixture.quad_part());
  const std::vector<stk::mesh::Bucket*>& buckets =
    bulk_data.buckets(stk::mesh::Face);
  std::vector<stk::mesh::Bucket*> live_buckets;
  get_buckets(selector, buckets, live_buckets);

  for (std::vector<stk::mesh::Bucket*>::iterator itr = live_buckets.begin();
       itr != live_buckets.end(); ++itr) {
    stk::mesh::Bucket& b = **itr;
    for (size_t entity_itr = 0; entity_itr < b.size(); ++entity_itr) {
      stk::mesh::Entity& entity = b[entity_itr];
      for (std::vector<unsigned>::const_iterator killed_itr = entity_ids_to_kill.begin(); killed_itr != entity_ids_to_kill.end(); ++killed_itr) {
        if (static_cast<unsigned>(entity.identifier()) == *killed_itr) {
          passed = 0;
        }
      }
    }
  }

  stk::all_reduce(pm, stk::ReduceMin<1>(&passed));



  return passed != 0;
}

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

bool element_death_use_case(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_size(pm) > 1) {
    return true;
  }

  GridFixture fixture(pm);
  stk::mesh::BulkData& bulk_data = fixture.bulk_data();
  stk::mesh::MetaData& meta_data = fixture.meta_data();

  // Arbitrarily chose faces 9, 10 to kill; see GridFixture.cpp
  std::vector<unsigned> entity_ids_to_kill;
  entity_ids_to_kill.push_back(9);
  entity_ids_to_kill.push_back(10);

  std::vector<stk::mesh::Entity*> entities_to_kill;
  if (stk::parallel_machine_rank(pm) == 0) {
    for (std::vector<unsigned>::const_iterator itr = entity_ids_to_kill.begin();
         itr != entity_ids_to_kill.end(); ++itr) {
      entities_to_kill.push_back(bulk_data.get_entity(stk::mesh::Face,
                                                      *itr));
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
    const stk::mesh::Entity* live_entity = itr->second.first;
    if (live_entity &&
        live_entity->owner_rank() == bulk_data.parallel_rank()) {
      ++new_side_count;
    }
  }

  // Ask for new globally unique side ids
  std::vector<size_t> requests(meta_data.entity_type_count(), 0);
  std::vector<stk::mesh::EntityKey> requested_keys;
  requests[stk::mesh::Edge] = new_side_count;
  bulk_data.generate_new_keys(requests, requested_keys);

  // Create boundaries between live and dead entities
  new_side_count = 0;
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
       itr != boundary.end(); ++itr) {
    stk::mesh::Entity* live_entity = itr->second.first;
    const unsigned local_id = itr->second.second;
    if (live_entity &&
        live_entity->owner_rank() == bulk_data.parallel_rank()) {
      stk::mesh::EntityId requested_id = 
        stk::mesh::entity_id(requested_keys[new_side_count]);
      stk::mesh::declare_element_side(bulk_data, requested_id, *live_entity,
                                      local_id,
                                      fixture.boundary_part());
      ++new_side_count;
    }
  }

  bulk_data.modification_end();

  // Check to see that results are correct
  // print live, print dead etc
  bool passed = true;
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
          passed = false;
        }
      }
    }
  }

  return passed;
}

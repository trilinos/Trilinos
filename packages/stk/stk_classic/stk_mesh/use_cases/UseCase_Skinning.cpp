/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <use_cases/UseCase_Skinning.hpp>

#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>

namespace {

static const size_t NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

void find_owned_nodes_with_relations_outside_closure(
    stk_classic::mesh::EntityVector & closure,
    stk_classic::mesh::Selector       select_owned,
    stk_classic::mesh::EntityVector & nodes)
{
  nodes.clear();

  //the closure is a sorted unique vector
  const stk_classic::mesh::EntityRank upward_rank = NODE_RANK + 1;
  const stk_classic::mesh::EntityId base_id = 0;
  stk_classic::mesh::EntityVector::iterator node_end = std::lower_bound(closure.begin(),
      closure.end(),
      stk_classic::mesh::EntityKey(upward_rank, base_id),
      stk_classic::mesh::EntityLess());

  for (stk_classic::mesh::EntityVector::iterator itr = closure.begin(); itr != node_end; ++itr) {
    stk_classic::mesh::Entity & node = **itr;

    if (select_owned(node)) {
      stk_classic::mesh::PairIterRelation relations_pair = node.relations();

      //loop over the relations and check to see if they are in the closure
      for (; relations_pair.first != relations_pair.second; ++relations_pair.first) {
        stk_classic::mesh::Entity * current_entity = (relations_pair.first->entity());

        //has relation outside of closure
        if ( !std::binary_search(node_end,
              closure.end(),
              current_entity,
              stk_classic::mesh::EntityLess()) )
        {
          nodes.push_back(&node);
          break;
        }
      }
    }
  }
}

void copy_nodes_and_break_relations( stk_classic::mesh::BulkData     & mesh,
                                     stk_classic::mesh::EntityVector & closure,
                                     stk_classic::mesh::EntityVector & nodes,
                                     stk_classic::mesh::EntityVector & new_nodes)
{
  for (size_t i = 0; i < nodes.size(); ++i) {
    stk_classic::mesh::Entity * entity = nodes[i];
    stk_classic::mesh::Entity * new_entity = new_nodes[i];

    stk_classic::mesh::PairIterRelation relations_pair = entity->relations();

    std::vector<stk_classic::mesh::EntitySideComponent> sides;

    //loop over the relations and check to see if they are in the closure
    for (; relations_pair.first != relations_pair.second;) {
      --relations_pair.second;
      stk_classic::mesh::Entity * current_entity = (relations_pair.second->entity());
      unsigned side_ordinal = relations_pair.second->identifier();

      if (stk_classic::mesh::in_receive_ghost(*current_entity)) {
        // TODO deleteing the ghost triggers a logic error at the
        // end of the NEXT modification cycle.  We need to fix this!
        //mesh.destroy_entity(current_entity);
        continue;
      }
      // check if has relation in closure
      else if ( std::binary_search(closure.begin(),
                                   closure.end(),
                                   current_entity,
                                   stk_classic::mesh::EntityLess()) )
      {
        sides.push_back(stk_classic::mesh::EntitySideComponent(current_entity,side_ordinal));
      }
    }

    //loop over the sides and break the relations between the old nodes
    //and set up the relations with the new
    for ( std::vector<stk_classic::mesh::EntitySideComponent>::iterator
          itr = sides.begin(); itr != sides.end(); ++itr) {
      mesh.destroy_relation(*(itr->entity), *entity, itr->side_ordinal);
      mesh.declare_relation(*(itr->entity), *new_entity, itr->side_ordinal);
    }

    //copy non-induced part membership from nodes[i] to new_nodes[i]
      //there are NO non-induced parts for this example

    //copy field data from nodes[i] to new_nodes[i]
    mesh.copy_entity_fields( *entity, *new_entity);

    if (entity->relations().empty()) {
      mesh.destroy_entity(entity);
    }

    if (new_entity->relations().empty()) {
      mesh.destroy_entity(new_entity);
    }


  }

}

void communicate_and_create_shared_nodes( stk_classic::mesh::BulkData & mesh,
                                          stk_classic::mesh::EntityVector   & nodes,
                                          stk_classic::mesh::EntityVector   & new_nodes)
{


  stk_classic::CommAll comm(mesh.parallel());

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk_classic::mesh::Entity & node = *nodes[i];
    stk_classic::mesh::Entity & new_node = *new_nodes[i];

    stk_classic::mesh::PairIterEntityComm entity_comm = node.sharing();

    for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

      unsigned proc = entity_comm.first->proc;
      comm.send_buffer(proc).pack<stk_classic::mesh::EntityKey>(node.key())
        .pack<stk_classic::mesh::EntityKey>(new_node.key());

    }
  }

  comm.allocate_buffers( mesh.parallel_size()/4 );

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk_classic::mesh::Entity & node = *nodes[i];
    stk_classic::mesh::Entity & new_node = *new_nodes[i];

    stk_classic::mesh::PairIterEntityComm entity_comm = node.sharing();

    for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

      unsigned proc = entity_comm.first->proc;
      comm.send_buffer(proc).pack<stk_classic::mesh::EntityKey>(node.key())
                            .pack<stk_classic::mesh::EntityKey>(new_node.key());

    }
  }

  comm.communicate();

  const stk_classic::mesh::PartVector no_parts;

  for (size_t process = 0; process < mesh.parallel_size(); ++process) {
    stk_classic::mesh::EntityKey old_key;
    stk_classic::mesh::EntityKey new_key;

    while ( comm.recv_buffer(process).remaining()) {

      comm.recv_buffer(process).unpack<stk_classic::mesh::EntityKey>(old_key)
                               .unpack<stk_classic::mesh::EntityKey>(new_key);

      stk_classic::mesh::Entity * old_entity = mesh.get_entity(old_key);
      stk_classic::mesh::Entity * new_entity = & mesh.declare_entity(new_key.rank(), new_key.id(), no_parts);

      nodes.push_back(old_entity);
      new_nodes.push_back(new_entity);

    }
  }
}

} // empty namespace

void separate_and_skin_mesh(
    stk_classic::mesh::fem::FEMMetaData & fem_meta,
    stk_classic::mesh::BulkData & mesh,
    stk_classic::mesh::Part     & skin_part,
    std::vector< stk_classic::mesh::EntityId > elements_to_separate,
    const stk_classic::mesh::EntityRank rank_of_element
    )
{
  stk_classic::mesh::EntityVector entities_to_separate;

  //select the entity only if the current process in the owner
  for (std::vector< stk_classic::mesh::EntityId>::const_iterator itr = elements_to_separate.begin();
      itr != elements_to_separate.end(); ++itr)
  {
    stk_classic::mesh::Entity * element = mesh.get_entity(rank_of_element, *itr);
    if (element != NULL && element->owner_rank() == mesh.parallel_rank()) {
      entities_to_separate.push_back(element);
    }
  }

  stk_classic::mesh::EntityVector entities_closure;
  stk_classic::mesh::find_closure(mesh,
      entities_to_separate,
      entities_closure);

  stk_classic::mesh::Selector select_owned = fem_meta.locally_owned_part();

  stk_classic::mesh::EntityVector nodes;
  find_owned_nodes_with_relations_outside_closure( entities_closure, select_owned, nodes);

  //ask for new nodes to represent the copies
  std::vector<size_t> requests(fem_meta.entity_rank_count(), 0);
  requests[NODE_RANK] = nodes.size();

  mesh.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk_classic::mesh::EntityVector new_nodes;
  mesh.generate_new_entities(requests, new_nodes);

  //communicate and create new nodes everywhere the old node is shared
  communicate_and_create_shared_nodes(mesh, nodes, new_nodes);

  copy_nodes_and_break_relations(mesh, entities_closure, nodes, new_nodes);

  mesh.modification_end();

  skin_mesh( mesh, rank_of_element, &skin_part);

  return;
}

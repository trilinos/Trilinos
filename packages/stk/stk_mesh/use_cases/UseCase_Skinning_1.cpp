/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <set>

#include <use_cases/HexFixture.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>

namespace {

  unsigned count_skin_entities( HexFixture & fixture) {

    const unsigned mesh_rank = stk::mesh::Element;

    stk::mesh::BulkData & mesh = fixture.m_bulk_data;
    stk::mesh::Part & skin_part = fixture.m_skin_part;

    stk::mesh::Selector select_skin = skin_part ;

    const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( mesh_rank -1);

    return count_selected_entities( select_skin, buckets);
  }

  void find_owned_nodes_with_relations_outside_closure(
      stk::mesh::EntityVector & closure,
      stk::mesh::Selector       select_owned,
      stk::mesh::EntityVector & nodes)
  {
    nodes.clear();

    //the closure is a sorted unique vector
    stk::mesh::EntityVector::iterator node_end = std::lower_bound(closure.begin(),
        closure.end(),
        stk::mesh::EntityKey(stk::mesh::Edge, 0),
        stk::mesh::EntityLess());

    for (stk::mesh::EntityVector::iterator itr = closure.begin(); itr != node_end; ++itr) {
      stk::mesh::Entity & node = **itr;

      if (select_owned(node)) {
        stk::mesh::PairIterRelation relations_pair = node.relations();

        //loop over the relations and check to see if they are in the closure
        for (; relations_pair.first != relations_pair.second; ++relations_pair.first) {
          stk::mesh::Entity * current_entity = (relations_pair.first->entity());

          //has relation outside of closure
          if ( !std::binary_search(node_end,
                closure.end(),
                current_entity,
                stk::mesh::EntityLess()) )
          {
            nodes.push_back(&node);
            break;
          }
        }
      }
    }
  }

  void copy_nodes_and_break_relations( stk::mesh::BulkData     & mesh,
                                       stk::mesh::EntityVector & closure,
                                       stk::mesh::EntityVector & nodes,
                                       stk::mesh::EntityVector & new_nodes)
  {


    for (size_t i = 0; i < nodes.size(); ++i) {
      stk::mesh::Entity & entity = *nodes[i];
      stk::mesh::PairIterRelation relations_pair = entity.relations();

      std::vector<stk::mesh::EntitySideComponent> sides;

      //loop over the relations and check to see if they are in the closure
      for (; relations_pair.first != relations_pair.second; ++relations_pair.first) {
        stk::mesh::Entity * current_entity = (relations_pair.first->entity());
        unsigned side_ordinal = relations_pair.first->identifier();

        //has relation in closure
        if ( std::binary_search(closure.begin(),
              closure.end(),
              current_entity,
              stk::mesh::EntityLess()) )
        {
          sides.push_back(stk::mesh::EntitySideComponent(current_entity,side_ordinal));
        }
      }

      //loop over the sides and break the relations between the old nodes
      //and set up the relations with the new
      for ( std::vector<stk::mesh::EntitySideComponent>::iterator itr = sides.begin();
           itr != sides.end(); ++itr)
      {
        mesh.destroy_relation(*(itr->entity), *nodes[i]);
        mesh.declare_relation(*(itr->entity), *new_nodes[i], itr->side_ordinal);
      }

      //copy non-induced part membership from nodes[i] to new_nodes[i]

      //copy field data from nodes[i] to new_nodes[i]

    }

  }

  void communicate_and_create_shared_nodes( stk::mesh::BulkData & mesh,
                                            stk::mesh::EntityVector   & nodes,
                                            stk::mesh::EntityVector   & new_nodes)
  {


    stk::CommAll comm(mesh.parallel());

    for (size_t i = 0; i < nodes.size(); ++i) {
      stk::mesh::Entity & node = *nodes[i];
      stk::mesh::Entity & new_node = *new_nodes[i];

      stk::mesh::PairIterEntityComm entity_comm = node.sharing();

      for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

        unsigned proc = entity_comm.first->proc;
        comm.send_buffer(proc).pack<stk::mesh::EntityKey>(node.key())
          .pack<stk::mesh::EntityKey>(new_node.key());

      }
    }

    comm.allocate_buffers( mesh.parallel_size()/4 );

    for (size_t i = 0; i < nodes.size(); ++i) {
      stk::mesh::Entity & node = *nodes[i];
      stk::mesh::Entity & new_node = *new_nodes[i];

      stk::mesh::PairIterEntityComm entity_comm = node.sharing();

      for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

        unsigned proc = entity_comm.first->proc;
        comm.send_buffer(proc).pack<stk::mesh::EntityKey>(node.key())
                              .pack<stk::mesh::EntityKey>(new_node.key());

      }
    }

    comm.communicate();

    const stk::mesh::PartVector no_parts;

    for (size_t process = 0; process < mesh.parallel_size(); ++process) {
      stk::mesh::EntityKey old_key;
      stk::mesh::EntityKey new_key;

      while ( comm.recv_buffer(process).remaining()) {

        comm.recv_buffer(process).unpack<stk::mesh::EntityKey>(old_key)
                                 .unpack<stk::mesh::EntityKey>(new_key);

        stk::mesh::Entity * old_entity = mesh.get_entity(old_key);
        stk::mesh::Entity * new_entity = & mesh.declare_entity(new_key.rank(), new_key.id(), no_parts);

        nodes.push_back(old_entity);
        new_nodes.push_back(new_entity);

      }
    }
  }
}

bool skinning_use_case_1(stk::ParallelMachine pm)
{
  const unsigned p_size = stk::parallel_machine_size( pm );
  if ( 1 < p_size ) { return true; }

  //setup the mesh
  HexFixture fixture(pm);

  stk::mesh::MetaData& meta_data = fixture.m_meta_data;
  stk::mesh::BulkData& mesh = fixture.m_bulk_data;

  stk::mesh::Part & skin_part = fixture.m_skin_part;

  stk::mesh::EntityVector entities_to_separate;

  stk::mesh::Entity * temp = fixture.m_elems[0][0][0];

  //select the entity only if the current process in the owner
  if (temp != NULL && temp->owner_rank() == mesh.parallel_rank()) {
    entities_to_separate.push_back(temp);
  }

  stk::mesh::EntityVector entities_closure;
  stk::mesh::find_closure(mesh,
      entities_to_separate,
      entities_closure);

  stk::mesh::Selector select_owned = meta_data.locally_owned_part();

  stk::mesh::EntityVector nodes;
  find_owned_nodes_with_relations_outside_closure( entities_closure, select_owned, nodes);

  //ask for new nodes to represent the copies
  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[stk::mesh::Node] = nodes.size();

  mesh.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  mesh.generate_new_entities(requests, new_nodes);

  //communicate and create new nodes everywhere the old node is shared
  communicate_and_create_shared_nodes(mesh, nodes, new_nodes);

  copy_nodes_and_break_relations(mesh, entities_closure, nodes, new_nodes);

  mesh.modification_end();

  skin_mesh( mesh, stk::mesh::Element, &skin_part);

  unsigned num_skin_entities = count_skin_entities(fixture);

  stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

  bool correct_skin = num_skin_entities == 60;

  return correct_skin;
}

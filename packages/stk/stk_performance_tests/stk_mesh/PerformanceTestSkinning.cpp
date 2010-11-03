/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FEMTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

namespace {

size_t count_skin_entities( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, stk::mesh::EntityRank skin_rank ) {

  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

void delete_skin( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, stk::mesh::EntityRank skin_rank ) {

  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  stk::mesh::EntityVector skin_entities;

  stk::mesh::get_selected_entities( select_skin, buckets, skin_entities);

  mesh.modification_begin();

  for ( stk::mesh::EntityVector::iterator i = skin_entities.begin(); i != skin_entities.end(); ++i) {
    mesh.destroy_entity(*i);
  }

  mesh.modification_end();


}

void update_skin( stk::mesh::BulkData & mesh, stk::mesh::Part *skin_part, stk::mesh::EntityRank element_rank ) {

  stk::mesh::EntityVector owned_elements, modified_elements;

  // select owned
  stk::mesh::Selector owned = mesh.mesh_meta_data().locally_owned_part();
  stk::mesh::get_selected_entities( owned,
                         mesh.buckets(element_rank),
                         owned_elements);

  for( stk::mesh::EntityVector::iterator i = owned_elements.begin();
      i != owned_elements.end(); ++i )
  {
    stk::mesh::Entity * entity= *i;
    if ( entity->log_query() == stk::mesh::EntityLogCreated ||
         entity->log_query() == stk::mesh::EntityLogModified ) {
     modified_elements.push_back(entity);
    }
  }


  stk::mesh::reskin_mesh(mesh, element_rank,
                         modified_elements, skin_part);
}

void find_owned_nodes_with_relations_outside_closure(
    stk::mesh::EntityVector & closure,
    stk::mesh::Selector       select_owned,
    stk::mesh::EntityVector & nodes)
{
  nodes.clear();

  //the closure is a sorted unique vector
  const stk::mesh::EntityRank upward_rank = stk::mesh::NodeRank + 1;
  stk::mesh::EntityVector::iterator node_end = std::lower_bound(closure.begin(),
      closure.end(),
      stk::mesh::EntityKey(upward_rank, 0),
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
    stk::mesh::Entity * entity = nodes[i];
    stk::mesh::Entity * new_entity = new_nodes[i];

    stk::mesh::PairIterRelation relations_pair = entity->relations();

    std::vector<stk::mesh::EntitySideComponent> sides;

    //loop over the relations and check to see if they are in the closure
    for (; relations_pair.first != relations_pair.second;) {
      --relations_pair.second;
      stk::mesh::Entity * current_entity = (relations_pair.second->entity());
      size_t side_ordinal = relations_pair.second->identifier();

      if (stk::mesh::in_receive_ghost(*current_entity)) {
        // TODO deleteing the ghost triggers a logic error at the
        // end of the NEXT modification cycle.  We need to fix this!
        //mesh.destroy_entity(current_entity);
        continue;
      }
      else if ( std::binary_search(closure.begin(),
            //has relation in closure
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
      mesh.destroy_relation(*(itr->entity), *entity);
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

      size_t proc = entity_comm.first->proc;
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

      size_t proc = entity_comm.first->proc;
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

void separate_and_skin_mesh(
    stk::mesh::MetaData & meta,
    stk::mesh::BulkData & mesh,
    stk::mesh::Part     & skin_part,
    stk::mesh::EntityVector & entities_to_separate
    )
{

  stk::mesh::EntityVector entities_closure;
  stk::mesh::find_closure(mesh,
      entities_to_separate,
      entities_closure);

  stk::mesh::Selector select_owned = meta.locally_owned_part();

  stk::mesh::EntityVector nodes;
  find_owned_nodes_with_relations_outside_closure( entities_closure, select_owned, nodes);

  //ask for new nodes to represent the copies
  std::vector<size_t> requests(meta.entity_rank_count(), 0);
  requests[stk::mesh::NodeRank] = nodes.size();

  mesh.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  mesh.generate_new_entities(requests, new_nodes);

  //communicate and create new nodes everywhere the old node is shared
  communicate_and_create_shared_nodes(mesh, nodes, new_nodes);

  copy_nodes_and_break_relations(mesh, entities_closure, nodes, new_nodes);

  mesh.modification_end();


  return;
}

}//end unnamped namespace

// \TODO Idea: ADD scaling test over mesh size and compute the slope.
// \TODO Idea: ADD different partitioning such that the reskinning spans more than one process.

STKUNIT_UNIT_TEST( PerformanceTestSkinning, large_cube)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const size_t p_size = stk::parallel_machine_size(pm);
  const size_t p_rank = stk::parallel_machine_rank(pm);

  // Every processor will be involved in detachment and skin-update up to 500 processors.
  // This puts 5000 elements on each process unless we are running with STL
  // in debug mode in which case we shrink the problem down in order
  // to keep things running in a reasonable amount of time.
#ifdef _GLIBCXX_DEBUG
  const size_t NX = p_size*10, NY = 4, NZ = 5;
#else
  const size_t NX = p_size*10, NY = 20, NZ = 25;
#endif

  /* timings[0] = create mesh
   * timings[1] = intial skin mesh
   * timings[2] = detach mesh
   * timings[3] = delete skin
   * timings[4] = reskin mesh
   * timings[5] = sum(timings[0:4])
   */
  double timings[6] = {0};
  double start_time = 0;

  //recreate skin
  for ( int test_run = 0; test_run < 4; ++test_run) {
    //create the mesh

    start_time = stk::wall_time();
    stk::mesh::fixtures::HexFixture fixture(pm,NX,NY,NZ);

    stk::mesh::MetaData & meta = fixture.m_meta_data;
    stk::mesh::BulkData & mesh = fixture.m_bulk_data;
    stk::mesh::TopologicalMetaData & top = fixture.m_top_data;

    stk::mesh::Part & skin_part = meta.declare_part("skin_part");
    meta.commit();

    fixture.generate_mesh();
    timings[0] = stk::wall_dtime(start_time);

    //intial skin of the mesh
    start_time = stk::wall_time();
    stk::mesh::skin_mesh(mesh, top.element_rank, &skin_part);
    timings[1] = stk::wall_dtime(start_time);

    stk::mesh::EntityVector entities_to_separate;

    if ( test_run < 2) {
      //detach 1/3 of the mesh
      size_t num_detached_this_proc = 0;
      for (size_t ix=NX/3; ix < 2*NX/3; ++ix) {
      for (size_t iy=0; iy < NY; ++iy) {
      for (size_t iz=0; iz < NZ; ++iz) {
        stk::mesh::Entity * element = fixture.elem(ix,iy,iz);
        if (element != NULL && element->owner_rank() == mesh.parallel_rank()) {
          entities_to_separate.push_back(element);
          num_detached_this_proc++;
        }
      }
      }
      }
      STKUNIT_EXPECT_TRUE( num_detached_this_proc > 0u );
    } else {
      //detach middle of the mesh
      for (size_t ix=NX/2; ix < NX/2+1; ++ix) {
      for (size_t iy=NY/2; iy < NY/2+1; ++iy) {
      for (size_t iz=NZ/2; iz < NZ/2+1; ++iz) {
        stk::mesh::Entity * element = fixture.elem(ix,iy,iz);
        if (element != NULL && element->owner_rank() == mesh.parallel_rank()) {
          entities_to_separate.push_back(element);
        }
      }
      }
      }
    }

    start_time = stk::wall_time();
    separate_and_skin_mesh(
        meta,
        mesh,
        skin_part,
        entities_to_separate
        );
    timings[2] = stk::wall_dtime(start_time);

    if (test_run%2 == 0) { // delete the skin
      start_time = stk::wall_time();
      delete_skin( mesh, skin_part, top.side_rank );
      timings[3] = stk::wall_dtime(start_time);

      //reskin the entire mesh
      start_time = stk::wall_time();
      stk::mesh::skin_mesh( mesh, top.element_rank, &skin_part);
      timings[4] = stk::wall_dtime(start_time);
    }
    else { //update the skin
      timings[3] = 0;

      //update the skin of the mesh
      start_time = stk::wall_time();
      update_skin( mesh, &skin_part, top.element_rank);
      timings[4] = stk::wall_dtime(start_time);
    }

    //total the timings
    timings[5] = 0;
    for (int i=0; i <5; ++i) {
      timings[5] += timings[i];
    }

    stk::all_reduce(pm, stk::ReduceMax<5>(timings));

    if (p_rank == 0) {
      std::cout << "\n\n";
      switch (test_run) {
        case 0:
          std::cout << "Recreate entire skin after detaching 1/3 of the mesh:\n";
          break;
        case 1:
          std::cout << "Update skin after detaching 1/3 of the mesh:\n";
          break;
        case 2:
          std::cout << "Recreate entire skin after detaching middle of the mesh:\n";
          break;
        case 3:
          std::cout << "Update skin after detaching middle of the mesh:\n";
          break;
      }

      std::cout << "Num procs: " << p_size << "\n";
      std::cout << "Mesh size: " << NX << 'x' << NY << 'x' << NZ << " = " << NX*NY*NZ << " elements\n";
      std::cout << "Total time: "     << timings[5] << "\n";
      std::cout << "\tCreate mesh: "  << timings[0] << "\n";
      std::cout << "\tInitial skin: " << timings[1] << "\n";
      std::cout << "\tDetach mesh: "  << timings[2] << "\n";
      std::cout << "\tDelete skin: "  << timings[3] << "\n";
      std::cout << "\tReskin:      "  << timings[4] << "\n";
      std::cout << "\n\n";
    }

    size_t num_skin_entities = count_skin_entities(mesh, skin_part, top.side_rank );

    stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

    size_t expected_num_skin;

    if ( test_run < 2) {
      expected_num_skin = 2*(NX*NY + NX*NZ + 3*NY*NZ);
    } else {
      expected_num_skin = 2*(NX*NY + NX*NZ + NY*NZ) + 12;
    }

    STKUNIT_EXPECT_EQUAL( num_skin_entities, expected_num_skin );
  }

}

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <stdexcept>

#include <gtest/gtest.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_unit_test_utils/timer.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_io/WriteMesh.hpp>

namespace {

using stk::mesh::EntityRank;

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

size_t count_skin_entities( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, EntityRank skin_rank ) {

  const stk::mesh::MetaData & fem_meta = mesh.mesh_meta_data();

  stk::mesh::Selector select_skin = skin_part & fem_meta.locally_owned_part()  ;

  const stk::mesh::BucketVector& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

void delete_skin( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, EntityRank skin_rank ) {

  const stk::mesh::MetaData & fem_meta = mesh.mesh_meta_data();

  stk::mesh::Selector select_skin = skin_part & fem_meta.locally_owned_part()  ;

  const stk::mesh::BucketVector& buckets = mesh.buckets( skin_rank );

  stk::mesh::EntityVector skin_entities;

  stk::mesh::get_selected_entities( select_skin, buckets, skin_entities);

  mesh.modification_begin();

  for ( stk::mesh::EntityVector::iterator i = skin_entities.begin(); i != skin_entities.end(); ++i) {
    mesh.destroy_entity(*i);
  }

  mesh.modification_end();
}

void update_skin( stk::mesh::BulkData & mesh, stk::mesh::Part *skin_part, EntityRank element_rank ) {

  stk::mesh::EntityVector owned_elements, modified_elements;

  // select owned
  const stk::mesh::MetaData & fem_meta = mesh.mesh_meta_data();
  stk::mesh::Selector owned = fem_meta.locally_owned_part();
  stk::mesh::get_selected_entities( owned,
                         mesh.buckets(element_rank),
                         owned_elements);

  for( stk::mesh::EntityVector::iterator i = owned_elements.begin();
      i != owned_elements.end(); ++i )
  {
    stk::mesh::Entity entity= *i;
    if ( mesh.state(entity) == stk::mesh::Created ||
         mesh.state(entity) == stk::mesh::Modified ) {
     modified_elements.push_back(entity);
    }
  }

  stk::mesh::PartVector add_parts(1, skin_part);
  stk::mesh::skin_mesh(mesh, add_parts);
}

void find_owned_nodes_with_relations_outside_closure(
    const stk::mesh::BulkData& mesh,
    stk::mesh::EntityVector & closure,
    stk::mesh::Selector       select_owned,
    stk::mesh::EntityVector & nodes)
{
  nodes.clear();

  //the closure is a sorted unique vector
  const EntityRank upward_rank = stk::topology::EDGE_RANK;
  stk::mesh::EntityVector::iterator node_end = std::lower_bound(closure.begin(),
      closure.end(),
      stk::mesh::EntityKey(upward_rank, 0),
      stk::mesh::EntityLess(mesh));

  for (stk::mesh::EntityVector::iterator itr = closure.begin(); itr != node_end; ++itr) {
    stk::mesh::Entity node = *itr;

    if (select_owned(mesh.bucket(node))) {

      const stk::mesh::Bucket &node_bkt = mesh.bucket(node);
      const stk::mesh::Ordinal node_bkt_idx = mesh.bucket_ordinal(node);

      //loop over the relations and check to see if they are in the closure
      for (EntityRank irank = stk::topology::BEGIN_RANK;
            irank != stk::topology::END_RANK;
            ++irank)
      {
        stk::mesh::Entity const * rel_entity_i   = node_bkt.begin(node_bkt_idx, irank);
        stk::mesh::Entity const * rel_entity_end = node_bkt.end(node_bkt_idx, irank);
        for ( ; rel_entity_i != rel_entity_end ; ++rel_entity_i ) {
          stk::mesh::Entity current_entity = *rel_entity_i;
          //has relation outside of closure
          if ( !std::binary_search(node_end,
              closure.end(),
              current_entity,
              stk::mesh::EntityLess(mesh)) )
          {
            nodes.push_back(node);
            break;
          }
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
    stk::mesh::Entity entity = nodes[i];
    stk::mesh::Entity new_entity = new_nodes[i];

    std::vector<stk::mesh::EntitySideComponent> sides;

    //loop over the relations and check to see if they are in the closure
    for (stk::mesh::EntityRank irank = stk::topology::END_RANK;
          irank != stk::topology::BEGIN_RANK; )
    {
      --irank;

      int num_rels = mesh.num_connectivity(entity, irank);
      stk::mesh::Entity const *rel_ents = mesh.begin(entity, irank);
      stk::mesh::ConnectivityOrdinal const *rel_ords = mesh.begin_ordinals(entity, irank);

      for (int k = num_rels - 1; k >= 0; --k)
      {
        stk::mesh::Entity current_entity = rel_ents[k];
        size_t side_ordinal = rel_ords[k];

        if (mesh.in_receive_ghost(mesh.entity_key(current_entity))) {
          // TODO deleteing the ghost triggers a logic error at the
          // end of the NEXT modification cycle.  We need to fix this!
          //mesh.destroy_entity(current_entity);
          continue;
        }
        else if ( std::binary_search(closure.begin(),
            //has relation in closure
            closure.end(),
            current_entity,
            stk::mesh::EntityLess(mesh)) )
        {
          sides.push_back(stk::mesh::EntitySideComponent(current_entity,side_ordinal));
        }
      }
    }

    //loop over the sides and break the relations between the old nodes
    //and set up the relations with the new
    for ( std::vector<stk::mesh::EntitySideComponent>::iterator itr = sides.begin();
        itr != sides.end(); ++itr)
    {
      mesh.destroy_relation(itr->entity, entity, itr->side_ordinal);
      mesh.declare_relation(itr->entity, new_entity, itr->side_ordinal);
    }

    //copy non-induced part membership from nodes[i] to new_nodes[i]
    //there are NO non-induced parts for this example

    //copy field data from nodes[i] to new_nodes[i]
    mesh.copy_entity_fields( entity, new_entity);

    if (mesh.has_no_relations(entity)) {
      mesh.destroy_entity(entity);
    }

    if (mesh.has_no_relations(new_entity)) {
      mesh.destroy_entity(new_entity);
    }
  }
}

void communicate_and_create_shared_nodes( stk::mesh::BulkData & mesh,
    stk::mesh::EntityVector   & nodes,
    stk::mesh::EntityVector   & new_nodes)
{


  stk::CommSparse comm(mesh.parallel());

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk::mesh::Entity node = nodes[i];
    stk::mesh::Entity new_node = new_nodes[i];

    std::vector<int> shared_procs;
    mesh.comm_shared_procs(mesh.entity_key(node),shared_procs);

    for (size_t shared_procs_i=0 ; shared_procs_i<shared_procs.size() ; ++shared_procs_i) {

      size_t proc = shared_procs[shared_procs_i];
      comm.send_buffer(proc).pack<stk::mesh::EntityKey>(mesh.entity_key(node))
        .pack<stk::mesh::EntityKey>(mesh.entity_key(new_node));

    }
  }

  comm.allocate_buffers();

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk::mesh::Entity node = nodes[i];
    stk::mesh::Entity new_node = new_nodes[i];

    std::vector<int> shared_procs;
    mesh.comm_shared_procs(mesh.entity_key(node),shared_procs);

    for (size_t shared_procs_i=0 ; shared_procs_i<shared_procs.size() ; ++shared_procs_i) {

      size_t proc = shared_procs[shared_procs_i];
      comm.send_buffer(proc).pack<stk::mesh::EntityKey>(mesh.entity_key(node))
        .pack<stk::mesh::EntityKey>(mesh.entity_key(new_node));

    }
  }

  comm.communicate();

  const stk::mesh::PartVector no_parts;

  for (int process = 0; process < mesh.parallel_size(); ++process) {
    stk::mesh::EntityKey old_key;
    stk::mesh::EntityKey new_key;

    while ( comm.recv_buffer(process).remaining()) {

      comm.recv_buffer(process).unpack<stk::mesh::EntityKey>(old_key)
        .unpack<stk::mesh::EntityKey>(new_key);

      stk::mesh::Entity old_entity = mesh.get_entity(old_key);

      EXPECT_EQ(stk::topology::NODE_RANK, new_key.rank());
      stk::mesh::Entity new_entity = mesh.declare_node(new_key.id(), no_parts);

      nodes.push_back(old_entity);
      new_nodes.push_back(new_entity);

    }
  }
}

void separate_and_skin_mesh(
    stk::mesh::MetaData & fem_meta,
    stk::mesh::BulkData & mesh,
    stk::mesh::Part     & skin_part,
    stk::mesh::EntityVector & entities_to_separate
    )
{

  stk::mesh::EntityVector entities_closure;
  stk::mesh::find_closure(mesh,
      entities_to_separate,
      entities_closure);

  stk::mesh::Selector select_owned = fem_meta.locally_owned_part();

  stk::mesh::EntityVector nodes;
  find_owned_nodes_with_relations_outside_closure( mesh, entities_closure, select_owned, nodes);

  //ask for new nodes to represent the copies
  std::vector<size_t> requests(fem_meta.entity_rank_count(), 0);
  requests[NODE_RANK] = nodes.size();

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

TEST( skinning_large_cube, skinning_large_cube)
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
  const size_t NX = p_size*100, NY = 100, NZ = 100;
#endif

  static const int TIMER_COUNT = 6;
  /* timings[0] = create mesh
   * timings[1] = intial skin mesh
   * timings[2] = detach mesh
   * timings[3] = delete skin
   * timings[4] = reskin mesh
   * timings[5] = sum(timings[0:4])
   */
  double timings[TIMER_COUNT] = {0};
  double timing_sums[TIMER_COUNT] = {0};
  const char* timer_names[TIMER_COUNT] = {"Create mesh", "Initial skin", "Detach mesh", "Delete skin", "Reskin", "Total time"};
  double start_time = 0;
  size_t memory_max[2] = {0};
  const char* memory_names[2] = {"Current memory", "Memory high water"};

  //recreate skin
  for ( int test_run = 0; test_run < 4; ++test_run) {
    //create the mesh

    start_time = stk::wall_time();
    stk::mesh::fixtures::HexFixture fixture(pm,NX,NY,NZ);
    const EntityRank element_rank = stk::topology::ELEMENT_RANK;
    const EntityRank side_rank = fixture.m_meta.side_rank();

    stk::mesh::MetaData & fem_meta = fixture.m_meta;
    stk::mesh::BulkData & mesh = fixture.m_bulk_data;

    stk::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();
    timings[0] = stk::wall_dtime(start_time);

    //intial skin of the mesh
    start_time = stk::wall_time();
    {
      stk::mesh::PartVector add_parts(1,&skin_part);
      stk::mesh::skin_mesh(mesh, add_parts);
    }
    timings[1] = stk::wall_dtime(start_time);

    stk::mesh::EntityVector entities_to_separate;

    if ( test_run < 2) {
      //detach 1/3 of the mesh
      size_t num_detached_this_proc = 0;
      for (size_t ix=NX/3; ix < 2*NX/3; ++ix) {
      for (size_t iy=0; iy < NY; ++iy) {
      for (size_t iz=0; iz < NZ; ++iz) {
        stk::mesh::Entity element = fixture.elem(ix,iy,iz);
        if (mesh.is_valid(element) && mesh.parallel_owner_rank(element) == mesh.parallel_rank()) {
          entities_to_separate.push_back(element);
          num_detached_this_proc++;
        }
      }
      }
      }
      EXPECT_TRUE( num_detached_this_proc > 0u );
    } else {
      //detach middle of the mesh
      for (size_t ix=NX/2; ix < NX/2+1; ++ix) {
      for (size_t iy=NY/2; iy < NY/2+1; ++iy) {
      for (size_t iz=NZ/2; iz < NZ/2+1; ++iz) {
        stk::mesh::Entity element = fixture.elem(ix,iy,iz);
        if (mesh.is_valid(element) && mesh.parallel_owner_rank(element) == mesh.parallel_rank()) {
          entities_to_separate.push_back(element);
        }
      }
      }
      }
    }

    start_time = stk::wall_time();
    separate_and_skin_mesh(
        fem_meta,
        mesh,
        skin_part,
        entities_to_separate
        );
    timings[2] = stk::wall_dtime(start_time);

    if (test_run%2 == 0) { // delete the skin
      start_time = stk::wall_time();
      delete_skin( mesh, skin_part, side_rank );
      timings[3] = stk::wall_dtime(start_time);

      //reskin the entire mesh
      start_time = stk::wall_time();
      {
        stk::mesh::PartVector add_parts(1,&skin_part);
        stk::mesh::skin_mesh( mesh, add_parts);
      }
      timings[4] = stk::wall_dtime(start_time);
    }
    else { //update the skin
      timings[3] = 0;

      //update the skin of the mesh
      start_time = stk::wall_time();
      update_skin( mesh, &skin_part, element_rank);
      timings[4] = stk::wall_dtime(start_time);
    }

    //total the timings
    timings[5] = 0;
    for (int i=0; i <5; ++i) {
      timings[5] += timings[i];
    }

    size_t mem_now = 0, mem_hwm = 0;
    stk::get_memory_usage(mem_now, mem_hwm);

    stk::all_reduce(pm, stk::ReduceMax<5>(timings));
    stk::all_reduce(pm, stk::ReduceMax<1>(&mem_now));
    stk::all_reduce(pm, stk::ReduceMax<1>(&mem_hwm));

    if (mem_now > memory_max[0]) {
      memory_max[0] = mem_now;
    }
    if (mem_hwm > memory_max[1]) {
      memory_max[1] = mem_hwm;
    }

    //compute sums
    for (int i=0; i<TIMER_COUNT; ++i) {
      timing_sums[i] += timings[i];
    }

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

    size_t num_skin_entities = count_skin_entities(mesh, skin_part, side_rank );

    stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

    size_t expected_num_skin = 0;

    if ( test_run < 2) {
      expected_num_skin = 2*(NX*NY + NX*NZ + 3*NY*NZ);
    } else {
      expected_num_skin = 2*(NX*NY + NX*NZ + NY*NZ) + 12;
    }

    EXPECT_EQ( num_skin_entities, expected_num_skin );
  }

  if (p_rank == 0) {
    stk::print_timers_and_memory(&timer_names[0], &timing_sums[0], TIMER_COUNT, &memory_names[0], &memory_max[0], 2);
  }

  stk::parallel_print_time_without_output_and_hwm(pm, timings[5]);
}

double run_skinning_large_cube_test(bool createEdges, unsigned numRuns, std::vector<size_t>& dims)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  EXPECT_TRUE(dims.size() >= 3);
  const size_t NX = dims[0], NY = dims[1], NZ = dims[2];

  double skinningTime = 0.0;

  for (unsigned testRun = 0; testRun < numRuns; ++testRun) {
    stk::mesh::fixtures::HexFixture fixture(pm,NX,NY,NZ);

    stk::mesh::MetaData & fem_meta = fixture.m_meta;

    stk::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();
    
    if(createEdges) {
      stk::mesh::create_edges(fixture.m_bulk_data);
    }

    stk::mesh::PartVector add_parts(1,&skin_part);
    double startTime = stk::wall_time();
    stk::mesh::skin_mesh(fixture.m_bulk_data, add_parts);
    skinningTime += stk::wall_dtime(startTime);
  }

  return skinningTime;
}

TEST(skinning_large_cube_perf_test, skinning_large_cube)
{
  std::vector<size_t> dims = {50, 50, 50};
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 5;
  stk::unit_test_util::BatchTimer batchTimer(MPI_COMM_WORLD);
  batchTimer.initialize_batch_timer();
  for (unsigned j = 0; j < NUM_RUNS; j++) {
    batchTimer.start_batch_timer();
    run_skinning_large_cube_test(true, NUM_ITERS, dims);
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);

}

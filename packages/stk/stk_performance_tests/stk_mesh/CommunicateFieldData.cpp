// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <gtest/gtest.h>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <vector>
#include <string>

#include <stk_unit_test_utils/getOption.h>

namespace
{

double initial_value1[3] = {-1, 2, -0.3};

void createNodalVectorField(stk::mesh::MetaData& meshMetaData, const std::string &field_name)
{
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &field1 = meshMetaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, field_name);
    stk::mesh::put_field_on_entire_mesh_with_initial_value(field1, initial_value1);
}

void createNodalVectorFields(stk::mesh::MetaData& meshMetaData)
{
    createNodalVectorField(meshMetaData, "disp");
    createNodalVectorField(meshMetaData, "vel");
    createNodalVectorField(meshMetaData, "acc");
    createNodalVectorField(meshMetaData, "force");
    meshMetaData.commit();
}

size_t setFieldData(stk::mesh::BulkData& mesh, stk::mesh::Selector select)
{
    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields(stk::topology::NODE_RANK);
    const stk::mesh::BucketVector& node_buckets = mesh.get_buckets(stk::topology::NODE_RANK, select);
    size_t num_nodes = 0;
    for(size_t i=0; i<node_buckets.size(); ++i) {
        stk::mesh::Bucket& bucket = *node_buckets[i];
        num_nodes += bucket.size();
        for(size_t f=0; f<fields.size(); ++f) {
            double* data = static_cast<double*>(stk::mesh::field_data(*fields[f], bucket));
            unsigned flen = stk::mesh::field_scalars_per_entity(*fields[f], bucket);
            for(size_t n=0; n<flen*bucket.size(); ++n) {
                data[n] = static_cast<double>(f);
            }
        }
    }

    return num_nodes;
}

void createMetaAndBulkData(stk::io::StkMeshIoBroker &exodusFileReader, std::string genMeshSpec = std::string("generated:100x100x100"))
{
    std::string exodusFileName = unitTestUtils::getOption("-i", "NO_FILE_SPECIFIED");
    if (exodusFileName == "NO_FILE_SPECIFIED") {
      exodusFileName = genMeshSpec;
    }

    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int myProc = stk::parallel_machine_rank(communicator);

    if (myProc==0) {
      std::cerr << "Starting To Read Mesh: " << exodusFileName << std::endl;
    }

    exodusFileReader.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
    createNodalVectorFields(stkMeshMetaData);

    Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data(new stk::mesh::BulkData(stkMeshMetaData, MPI_COMM_WORLD, false, NULL));
    exodusFileReader.set_bulk_data(arg_bulk_data);
    stk::mesh::BulkData& stkMeshBulkData = *arg_bulk_data;

    bool delay_field_data_allocation = true;
    exodusFileReader.populate_mesh(delay_field_data_allocation);
    exodusFileReader.populate_field_data();

    size_t num_ghostings = stkMeshBulkData.ghostings().size();
    std::vector<size_t> num_ghost_nodes(num_ghostings,0);
    for (size_t ghost_i = 0 ; ghost_i < num_ghostings ; ++ghost_i) {
        num_ghost_nodes[ghost_i] = setFieldData(stkMeshBulkData, stkMeshBulkData.ghosting_part(*stkMeshBulkData.ghostings()[ghost_i]));
    }

    for(int i=0; i<stkMeshBulkData.parallel_size(); ++i) {
        stk::parallel_machine_barrier(MPI_COMM_WORLD);
        if (stkMeshBulkData.parallel_rank() == i) {
            for (size_t ghost_i = 0 ; ghost_i<num_ghostings ; ++ghost_i) {
                std::cerr <<"proc "<<stkMeshBulkData.parallel_rank()<<", Number of " << stkMeshBulkData.ghostings()[ghost_i]->name() << " nodes: " << num_ghost_nodes[ghost_i] << std::endl;
            }
        }

        stk::parallel_machine_barrier(MPI_COMM_WORLD);
        std::cerr << std::flush;
        stk::parallel_machine_barrier(MPI_COMM_WORLD);
    }
    if (myProc==0) {
      std::cerr << "Finished Reading Mesh" << std::endl;
    }

    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    size_t numElements = entityCounts[stk::topology::ELEMENT_RANK];
    size_t numNodes = entityCounts[stk::topology::NODE_RANK];

    if (myProc==0) {
      std::cerr << "Number of elements: " << numElements << std::endl;
      std::cerr << "Number of nodes: " << numNodes << std::endl;
    }
}

void test_communicate_field_data_all_ghosting(stk::mesh::BulkData& mesh, int num_iters)
{
    const int my_proc = mesh.parallel_rank();
    if (my_proc == 0) {
        std::cerr << "Calling communicate_field_data " << num_iters << " times"<<std::endl;
    }

    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields();
    std::vector<const stk::mesh::FieldBase*> const_fields(fields.size());
    for(size_t i=0; i<fields.size(); ++i) {
        const_fields[i] = fields[i];
    }

    double start_time = stk::cpu_time();

    for(int iter=0; iter<num_iters; ++iter) {
        for (size_t ghost_i = 0 ; ghost_i<mesh.ghostings().size() ; ++ghost_i) {
            stk::mesh::communicate_field_data(*mesh.ghostings()[ghost_i], const_fields);
        }
    }

    double stk_comm_time_with_loop = stk::cpu_time() - start_time;

    start_time = stk::cpu_time();
    for(int iter=0; iter<num_iters; ++iter) {
        stk::mesh::communicate_field_data(mesh, const_fields);
    }

    double stk_comm_time_bulk = stk::cpu_time() - start_time;

    double max_time_with_loop=0;
    double max_time_bulk=0;
    MPI_Reduce(static_cast<void*>(&stk_comm_time_with_loop), static_cast<void*>(&max_time_with_loop), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, mesh.parallel());
    MPI_Reduce(static_cast<void*>(&stk_comm_time_bulk), static_cast<void*>(&max_time_bulk), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, mesh.parallel());

    if ( my_proc == 0 ) {
      std::cerr << "Time to do communicate_field_data with loop: " << max_time_with_loop << std::endl;
      std::cerr << "Time to do communicate_field_data with bulk: " << max_time_bulk << std::endl;
    }

    stk::parallel_print_time_for_performance_compare(mesh.parallel(), stk_comm_time_bulk);
}

void test_communicate_field_data_ghosting(stk::mesh::BulkData& mesh, const stk::mesh::Ghosting& ghosting, int num_iters)
{
    const int my_proc = mesh.parallel_rank();
    if (my_proc == 0) {
        std::cerr << "Calling communicate_field_data " << num_iters << " times"<<std::endl;
    }

    const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
    const stk::mesh::FieldVector& fields = meta.get_fields();
    std::vector<const stk::mesh::FieldBase*> const_fields(fields.size());
    for(size_t i=0; i<fields.size(); ++i) {
        const_fields[i] = fields[i];
    }

    double start_time = stk::cpu_time();

    for(int iter=0; iter<num_iters; ++iter) {
        stk::mesh::communicate_field_data(ghosting, const_fields);
    }

    double stk_comm_time = stk::cpu_time() - start_time;

    double max_time=0;
    MPI_Reduce(static_cast<void*>(&stk_comm_time), static_cast<void*>(&max_time), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, MPI_COMM_WORLD);

    if ( my_proc == 0 ) {
      std::cerr << "Time to do communicate_field_data: " << max_time << std::endl;
    }
    stk::parallel_print_time_for_performance_compare(mesh.parallel(), stk_comm_time);
}

void addPartToGhosting(stk::mesh::BulkData & bulk, const std::string & partName, stk::mesh::Ghosting& ghost)
{
    const int numProcs = bulk.parallel_size();
    const int myProc = bulk.parallel_rank();
    stk::mesh::MetaData &meta = bulk.mesh_meta_data();
    std::vector<stk::mesh::EntityProc> entities_to_ghost;
    stk::mesh::Selector surface_selector = *meta.get_part(partName) & meta.locally_owned_part();
    stk::mesh::EntityVector surface_nodes;
    stk::mesh::get_selected_entities(surface_selector,bulk.buckets(stk::topology::NODE_RANK), surface_nodes);
    const int next_proc = (myProc+1)%numProcs;
    for (size_t i=0 ; i<surface_nodes.size() ; ++i) {
        entities_to_ghost.push_back(stk::mesh::EntityProc(surface_nodes[i],next_proc));
    }
    bulk.change_ghosting( ghost, entities_to_ghost );
}

TEST(CommunicateFieldData, copy_to_all)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs < 2) {
      return;
    }

    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    std::string genMeshSpec = "generated:100x100x100|sideset:xXyY";
    createMetaAndBulkData(exodusFileReader,genMeshSpec);
    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();

    stkMeshBulkData.modification_begin();
    stk::mesh::Ghosting & ghosting_1 = stkMeshBulkData.create_ghosting( "CUSTOM_1" );
    stk::mesh::Ghosting & ghosting_2 = stkMeshBulkData.create_ghosting( "CUSTOM_2" );
    stk::mesh::Ghosting & ghosting_3 = stkMeshBulkData.create_ghosting( "CUSTOM_3" );
    stk::mesh::Ghosting & ghosting_4 = stkMeshBulkData.create_ghosting( "CUSTOM_4" );

    addPartToGhosting(stkMeshBulkData, "surface_1", ghosting_1);
    addPartToGhosting(stkMeshBulkData, "surface_2", ghosting_2);
    addPartToGhosting(stkMeshBulkData, "surface_3", ghosting_3);
    addPartToGhosting(stkMeshBulkData, "surface_4", ghosting_4);
    stkMeshBulkData.modification_end();


    size_t num_ghostings = stkMeshBulkData.ghostings().size();
    std::ostringstream oss;
    for (size_t ghost_i = 0 ; ghost_i < num_ghostings ; ++ghost_i) {
        size_t num_nodes = stk::mesh::count_selected_entities(stkMeshBulkData.ghosting_part(*stkMeshBulkData.ghostings()[ghost_i]),stkMeshBulkData.buckets(stk::topology::NODE_RANK));
        oss <<"proc "<<stkMeshBulkData.parallel_rank()<<", Number of " << stkMeshBulkData.ghostings()[ghost_i]->name() << " nodes: " << num_nodes << std::endl;
    }

    std::cerr << oss.str() << std::endl;

    test_communicate_field_data_all_ghosting(stkMeshBulkData, 1000);
}


TEST(CommunicateFieldData, Ghosting)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs < 2) {
      return;
    }
  
    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    createMetaAndBulkData(exodusFileReader);

    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();
    const stk::mesh::Ghosting& aura_ghosting = stkMeshBulkData.aura_ghosting();
    test_communicate_field_data_ghosting(stkMeshBulkData, aura_ghosting, 1000);
}

}

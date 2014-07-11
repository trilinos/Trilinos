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
#include <stk_util/parallel/Parallel.hpp>

#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <vector>
#include <string>

#include <optionParsing/getOption.h>

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
    const stk::mesh::FieldVector& fields = meta.get_fields();
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

void createMetaAndBulkData(stk::io::StkMeshIoBroker &exodusFileReader)
{
    std::string exodusFileName = unitTestUtils::getOption("-i", "NO_FILE_SPECIFIED");
    if (exodusFileName == "NO_FILE_SPECIFIED") {
      exodusFileName = "generated:5x5x5";
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

    size_t num_shared_nodes = setFieldData(stkMeshBulkData, stkMeshMetaData.globally_shared_part());
    size_t num_aura_nodes = setFieldData(stkMeshBulkData, stkMeshMetaData.aura_part());

    for(int i=0; i<stkMeshBulkData.parallel_size(); ++i) {
        stk::parallel_machine_barrier(MPI_COMM_WORLD);
        if (stkMeshBulkData.parallel_rank() == i) {
            std::cerr <<"proc "<<stkMeshBulkData.parallel_rank()<<", Number of shared nodes: " << num_shared_nodes << std::endl;
            std::cerr <<"proc "<<stkMeshBulkData.parallel_rank()<<", Number of aura nodes: " << num_aura_nodes << std::endl;
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

    double max_time;
    MPI_Reduce(static_cast<void*>(&stk_comm_time), static_cast<void*>(&max_time), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, MPI_COMM_WORLD);

    if ( my_proc == 0 ) {
      std::cerr << "Time to do communicate_field_data: " << max_time << std::endl;
    }
}

void test_communicate_field_data_shared(stk::mesh::BulkData& mesh, int num_iters)
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
        stk::CommAll commall;
        stk::mesh::communicate_field_data(mesh, const_fields.size(), &const_fields[0], commall);
    }

    double stk_comm_time = stk::cpu_time() - start_time;

    double max_time;
    MPI_Reduce(static_cast<void*>(&stk_comm_time), static_cast<void*>(&max_time), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, MPI_COMM_WORLD);

    if ( my_proc == 0 ) {
      std::cerr << "Time to do communicate_field_data: " << max_time << std::endl;
    }
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
    const stk::mesh::Ghosting& aura_ghosting = *stkMeshBulkData.ghostings()[1];
    test_communicate_field_data_ghosting(stkMeshBulkData, aura_ghosting, 100);
}

TEST(CommunicateFieldData, Shared)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if (numProcs < 2) {
      return;
    }

    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    createMetaAndBulkData(exodusFileReader);

    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();

    test_communicate_field_data_shared(stkMeshBulkData, 100);
}

}

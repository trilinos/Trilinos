#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

namespace
{
TEST(StkMesh, UsingStkIO)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string exodusFileName = "example.exo";

    // Use STK IO to populate a STK Mesh
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.open_mesh_database(exodusFileName, stk::io::READ_MESH);
    // Create MetaData in STK IO
    exodusFileReader.create_input_mesh();
    // Create BulkData in STK IO
    exodusFileReader.populate_bulk_data();

    // Get the STK Mesh (BulkData and MetaData) from STK IO
    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();

    // Test if the STK Mesh has 256 elements. Other examples will discuss details below.
    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<stk::mesh::EntityRank> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(256u, entityCounts[stk::topology::ELEMENT_RANK]);
}
}

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
namespace stk { namespace mesh { class BulkData; } }

namespace
{
TEST(StkMesh, UsingStkIO)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string exodusFileName = "example.exo";

    // Use STK IO to populate a STK Mesh
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    // Create MetaData in STK IO
    exodusFileReader.create_input_mesh();
    // Create BulkData in STK IO
    exodusFileReader.populate_bulk_data();

    // Get the STK Mesh (BulkData and MetaData) from STK IO
    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();

    // Test if the STK Mesh has 256 elements. Other examples will discuss details below.
    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(256u, entityCounts[stk::topology::ELEMENT_RANK]);
}
}

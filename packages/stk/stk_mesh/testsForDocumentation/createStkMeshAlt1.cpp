
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH

namespace
{
TEST(StkMesh, CreateStkMesh)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    // Creation of STK Mesh objects.
    // MetaData creates the universal_part, locally-owned part, and globally shared part.
    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, communicator);

    // STK IO module will be described in separate chapter.
    // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
    // The order of the following lines in {} are important
    {
        stk::io::StkMeshIoBroker exodusFileReader(communicator);

        // Inform STK IO which STK Mesh objects to populate later
        exodusFileReader.set_bulk_data(stkMeshBulkData);

        const std::string exodusFileName = "example.exo";
        exodusFileReader.open_mesh_database(exodusFileName, stk::io::READ_MESH);

        // Populate the MetaData which has the descriptions of the Parts and Fields.
        exodusFileReader.create_input_mesh();

        // Populate entities in STK Mesh from Exodus file
        exodusFileReader.populate_bulk_data();
    }

    // Test if the STK Mesh has 256 elements. Other examples will discuss details below.
    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(256u, entityCounts[stk::topology::ELEMENT_RANK]);
}

}

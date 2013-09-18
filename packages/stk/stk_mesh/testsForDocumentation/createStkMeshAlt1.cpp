#include <gtest/gtest.h>

#include <string>
#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

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
        stk::io::MeshData exodusFileReader(communicator);

        // Inform STK IO which STK Mesh objects to populate later
        exodusFileReader.set_bulk_data(stkMeshBulkData);

        const std::string exodusFileName = "example.exo";
        exodusFileReader.open_mesh_database(exodusFileName);

        // Populate the MetaData which has the descriptions of the Parts and Fields.
        exodusFileReader.create_input_mesh();

        // Populate entities in STK Mesh from Exodus file
        exodusFileReader.populate_bulk_data();
    }

    // Test if the STK Mesh has 256 elements. Other examples will discuss details below.
    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<stk::mesh::EntityRank> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(256u, entityCounts[stk::topology::ELEMENT_RANK]);
}

}

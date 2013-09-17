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

const int spatialDim = 3;

TEST(StkMesh, CreateStkMesh)
{
    MPI_Comm communicator = MPI_COMM_NULL;
    stk::io::MeshData exodusFileReader(communicator);

    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, communicator);
    exodusFileReader.set_bulk_data(stkMeshBulkData);

    const std::string exodusFileName = "example.exo";
    exodusFileReader.open_mesh_database(exodusFileName);

    exodusFileReader.create_input_mesh();

    exodusFileReader.populate_bulk_data();

    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<stk::mesh::EntityRank> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(256u, entityCounts[stk::topology::ELEMENT_RANK]);
}

}

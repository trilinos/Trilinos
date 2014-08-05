#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <unistd.h>                     // for unlink
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
//-BEGIN    
  TEST(StkMeshHowTo, UsingStkIO)
  {
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string fileName = "generated:8x8x8";
    //The 'generated:' syntax in fileName causes a hex mesh to be generated in memory.
    //If an exodus file-name is used instead, then the mesh is read from file. The code
    //below remains the same in either case.

    // Use STK IO to populate a STK Mesh
    stk::io::StkMeshIoBroker meshReader(communicator);
    meshReader.add_mesh_database(fileName, stk::io::READ_MESH);
    // Create MetaData in STK IO
    meshReader.create_input_mesh();
    // Create BulkData in STK IO
    meshReader.populate_bulk_data();

    // Get the STK Mesh (BulkData and MetaData) from STK IO
    stk::mesh::MetaData &stkMeshMetaData = meshReader.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = meshReader.bulk_data();

    // Test if the STK Mesh has 512 elements. Other examples will discuss details below.
    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    EXPECT_EQ(512u, entityCounts[stk::topology::ELEMENT_RANK]);
    unlink(fileName.c_str());
  }
//-END
}

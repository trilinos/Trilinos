#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp> // for create_faces
#include <stk_mesh/base/CreateEdges.hpp> // for create_edges
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

  TEST(StkMeshHowTo, CreateFacesEdgesHex)
  {
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);

    const std::string generatedFileName = "generated:8x8x8";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    // ============================================================
    //+ EXAMPLE
    //+ Create the faces..
    stk::mesh::create_faces(stkIo.bulk_data());

    //+ Create the edges..
    stk::mesh::create_edges(stkIo.bulk_data());

    // ==================================================
    // VERIFICATION
    stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
    EXPECT_EQ( 512u, entityCounts[stk::topology::ELEMENT_RANK]);
    EXPECT_EQ(1728u, entityCounts[stk::topology::FACE_RANK]);
    EXPECT_EQ(1944u, entityCounts[stk::topology::EDGE_RANK]);
  }
}

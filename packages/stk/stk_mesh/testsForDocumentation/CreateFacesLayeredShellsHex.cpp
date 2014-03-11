#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp> // for create_faces
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

  TEST(StkMeshHowTo, CreateFacesLayeredShellsHex)
  {
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);

    // Generate a mesh containing 1 hex part and 12 shell parts
    // Shells are layered 2 deep.
    const std::string generatedFileName = "generated:8x8x8|shell:xxyyzzXYZXYZ";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    // ============================================================
    //+ EXAMPLE
    //+ Create the faces 
    stk::mesh::create_faces(stkIo.bulk_data());

    // ==================================================
    // VERIFICATION
    stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
    EXPECT_EQ(1280u, entityCounts[stk::topology::ELEMENT_RANK]);
    //+ The shell faces are the same as the boundary hex faces
    EXPECT_EQ(1728u, entityCounts[stk::topology::FACE_RANK]);

    // Edges are not generated, only faces.
    EXPECT_EQ(0u,   entityCounts[stk::topology::EDGE_RANK]);
  }
}

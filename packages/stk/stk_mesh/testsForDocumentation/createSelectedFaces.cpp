#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector
namespace stk { namespace mesh { class BulkData; } }

namespace
{

  TEST(StkMeshHowTo, CreateSelectedFacesHex)
  {
    // ============================================================
    // INITIALIZATION
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkIo(communicator);

    // Generate a mesh containing 1 hex part and 6 shell parts
    const std::string generatedFileName = "generated:8x8x8|shell:xyzXYZ";
    stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();
    const stk::mesh::PartVector &all_parts = stkIo.meta_data().get_mesh_parts();

    // ============================================================
    //+ EXAMPLE
    //+ Create a selector containing just the shell parts.
    stk::mesh::Selector shell_subset;
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];
      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
	shell_subset |= *part;
      }
    }

    //+ Create the faces on just the selected shell parts.
    stk::mesh::create_faces(stkIo.bulk_data(), shell_subset);

    // ==================================================
    // VERIFICATION
    stk::mesh::Selector allEntities = stkIo.meta_data().universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkIo.bulk_data(), entityCounts);
    EXPECT_EQ( 896u, entityCounts[stk::topology::ELEMENT_RANK]);
    EXPECT_EQ( 384u, entityCounts[stk::topology::FACE_RANK]);

    // Edges are not generated, only faces.
    EXPECT_EQ(0u,   entityCounts[stk::topology::EDGE_RANK]);
  }
}

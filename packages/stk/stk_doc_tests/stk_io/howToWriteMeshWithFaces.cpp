#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace
{

TEST(StkIoHowTo, WriteMeshWithFaceBlock)
{
  std::string filename = "output.exo";
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Part* part = &meta.declare_part_with_topology("faceBlock", stk::topology::QUAD_4);
    stk::io::put_face_block_io_part_attribute(*part);
    stk::io::fill_mesh("generated:1x1x1", *bulk);
    bool connectFacesToEdges = true;
    stk::mesh::create_all_sides(*bulk, meta.universal_part(), {part}, connectFacesToEdges);

    unsigned numFaces = stk::mesh::count_selected_entities(meta.universal_part(), bulk->buckets(stk::topology::FACE_RANK));
    EXPECT_EQ(6u, numFaces);

    stk::io::StkMeshIoBroker ioBroker;
    ioBroker.set_bulk_data(*bulk);
    stk::io::write_mesh(filename, ioBroker, stk::io::WRITE_RESULTS);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::io::fill_mesh(filename, *bulk);

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(*bulk, entityCounts);
    EXPECT_EQ(6u, entityCounts[stk::topology::FACE_RANK]);
    EXPECT_EQ(1u, entityCounts[stk::topology::ELEM_RANK]);

    stk::mesh::EntityVector faces;
    stk::mesh::get_entities(*bulk, stk::topology::EDGE_RANK, faces);
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(*bulk, stk::topology::ELEM_RANK, elems);

    for(stk::mesh::Entity face : faces) {
      const stk::mesh::Entity* faceElement = bulk->begin_elements(face);
      EXPECT_EQ(1u, bulk->num_elements(face));
      EXPECT_EQ(faceElement[0], elems[0]);
    }

    EXPECT_EQ(6u, bulk->num_faces(elems[0]));
  }

  unlink(filename.c_str());
}

}

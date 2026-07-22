#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>

namespace
{

TEST(StkIoHowTo, WriteMeshWithEdges)
{
  std::string filename = "output.exo";
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Part* edgeBlockPart = &meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);
    stk::io::put_edge_block_io_part_attribute(*edgeBlockPart);
    stk::io::fill_mesh("generated:1x1x1", *bulk);
    stk::mesh::create_edges(*bulk, meta.universal_part(), edgeBlockPart);

    unsigned numEdges = stk::mesh::count_selected_entities(meta.universal_part(), bulk->buckets(stk::topology::EDGE_RANK));
    EXPECT_EQ(12u, numEdges);

    stk::io::StkMeshIoBroker ioBroker;
    ioBroker.set_bulk_data(bulk);
    stk::io::write_mesh(filename, ioBroker, stk::io::WRITE_RESULTS);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();
    stk::io::fill_mesh(filename, *bulk);

    const stk::mesh::Part* edgeBlockPart = meta.get_part("edgeBlock");
    ASSERT_FALSE(edgeBlockPart == nullptr);
    EXPECT_TRUE(stk::io::is_part_io_part(*edgeBlockPart));
    EXPECT_TRUE(stk::io::is_part_edge_block_io_part(*edgeBlockPart));

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(*bulk, entityCounts);
    EXPECT_EQ(12u, entityCounts[stk::topology::EDGE_RANK]);
    EXPECT_EQ(1u, entityCounts[stk::topology::ELEM_RANK]);

    stk::mesh::EntityVector edges;
    stk::mesh::get_entities(*bulk, stk::topology::EDGE_RANK, edges);
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(*bulk, stk::topology::ELEM_RANK, elems);

    for(stk::mesh::Entity edge : edges) {
      const stk::mesh::Entity* edgeElement = bulk->begin_elements(edge);
      EXPECT_EQ(1u, bulk->num_elements(edge));
      EXPECT_EQ(edgeElement[0], elems[0]);
    }

    EXPECT_EQ(12u, bulk->num_edges(elems[0]));
  }

  unlink(filename.c_str());
}

TEST(StkIoHowTo, Write2DMeshWithEdges)
{
  std::string filename = "output2D.exo";
  {
    unsigned nx = 2, ny = 2;
    stk::mesh::fixtures::QuadFixture quadFixture(MPI_COMM_WORLD, nx, ny);
    stk::mesh::Part& sidesetPart = quadFixture.m_meta.declare_part("surface_1", stk::topology::EDGE_RANK);
    stk::mesh::Part& edgeBlockPart = quadFixture.m_meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);

    stk::io::put_io_part_attribute(sidesetPart);
    stk::io::put_edge_block_io_part_attribute(edgeBlockPart);
    stk::io::put_io_part_attribute(quadFixture.m_quad_part);

    quadFixture.generate_mesh();

    stk::mesh::BulkData& bulk = quadFixture.m_bulk_data;

    stk::mesh::create_edges(bulk, quadFixture.m_quad_part, &edgeBlockPart);
    stk::mesh::create_exposed_block_boundary_sides(bulk, quadFixture.m_quad_part, {&sidesetPart});

    unsigned numSidesetEdges = stk::mesh::count_selected_entities(sidesetPart, bulk.buckets(stk::topology::EDGE_RANK));
    unsigned expectedNumSidesetEdges = (nx + ny) * 2;
    EXPECT_EQ(expectedNumSidesetEdges, numSidesetEdges);

    unsigned numEdgeBlockEdges = stk::mesh::count_selected_entities(edgeBlockPart, bulk.buckets(stk::topology::EDGE_RANK));
    unsigned expectedNumEdgeBlockEdges = ny*(nx+1) + nx*(ny+1);
    EXPECT_EQ(expectedNumEdgeBlockEdges, numEdgeBlockEdges);

    stk::io::StkMeshIoBroker ioBroker;
    ioBroker.set_bulk_data(bulk);
    stk::io::write_mesh(filename, ioBroker, stk::io::WRITE_RESULTS);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();
    stk::io::fill_mesh(filename, *bulk);

    EXPECT_EQ(2u, meta.spatial_dimension());

    const stk::mesh::Part* edgeBlockPart = meta.get_part("edgeBlock");
    ASSERT_FALSE(edgeBlockPart == nullptr);
    EXPECT_TRUE(stk::io::is_part_io_part(*edgeBlockPart));
    EXPECT_TRUE(stk::io::is_part_edge_block_io_part(*edgeBlockPart));

    const stk::mesh::Part* sidesetPart = meta.get_part("surface_1");
    ASSERT_FALSE(sidesetPart == nullptr);
    EXPECT_TRUE(stk::io::is_part_io_part(*sidesetPart));
    EXPECT_FALSE(stk::io::is_part_edge_block_io_part(*sidesetPart));

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(*bulk, entityCounts);
    EXPECT_EQ(12u, entityCounts[stk::topology::EDGE_RANK]);
    EXPECT_EQ(4u, entityCounts[stk::topology::ELEM_RANK]);
  }

  unlink(filename.c_str());
}

}

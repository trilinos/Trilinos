#include <Akri_AllReduce.hpp>
#include <gtest/gtest.h>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_FacetMesh.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_ScopedTestFile.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_GroupingEntity.h>
#include <Ioss_IOFactory.h>
#include <Ioss_SubSystem.h>
#include <stk_io/IossBridge.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/EnvData.hpp>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

void put_node_with_id_into_nodeset(stk::mesh::BulkData & mesh, stk::mesh::Part & nodeset, const stk::mesh::EntityId id)
{
  const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, id);
  mesh.modification_begin();
  if (mesh.is_valid(node) && mesh.bucket(node).owned())
    mesh.change_entity_parts(node, stk::mesh::PartVector{&nodeset});
  mesh.modification_end();
}

void change_block_for_element_with_id(stk::mesh::BulkData & mesh, const stk::mesh::Part & oldBlock, const stk::mesh::Part & newBlock, const stk::mesh::EntityId id)
{
  const stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, id);
  mesh.modification_begin();
  if (mesh.is_valid(elem) && mesh.bucket(elem).owned())
    mesh.change_entity_parts(elem, stk::mesh::ConstPartVector{&newBlock}, stk::mesh::ConstPartVector{&oldBlock});
  mesh.modification_end();
}

int get_num_nodes_in_exodus_file(const std::string & filename)
{
  int numNodes = 0;
  if ( 0 == stk::EnvData::parallel_rank() )
  {
    /* open file */
    Ioss::DatabaseIO *db = Ioss::IOFactory::create("exodusII", filename.c_str(), Ioss::READ_MODEL, MPI_COMM_SELF);
    if ( !db ) {
      STK_ThrowErrorMsg("error reading file " << filename);
    }
    std::unique_ptr<Ioss::Region> io = std::make_unique<Ioss::Region>(db, "EXOSurface IC Region");

    numNodes = io->get_property("node_count").get_int();
  }
  krino::all_reduce_sum(stk::EnvData::parallel_comm(), numNodes);
  return numNodes;
}

void output_mesh_and_test_num_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const int goldNumOutputNodes)
{
  const krino::ScopedTestFile outputFile("outputTest.e");
  krino::output_composed_mesh_with_fields(mesh, outputSelector, outputFile.name(), 1, 0.);
  EXPECT_EQ(goldNumOutputNodes, get_num_nodes_in_exodus_file(outputFile.name()));
}

TEST(OutputUtils, createMesh_outputMeshWithSelectingBlock_allNodesOutput)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (krino::is_parallel_io_enabled() || stk::parallel_machine_size(comm) == 1)
  {
    krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, comm);
    krino::populate_bounding_box_mesh_and_activate(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 0.3333);

    stk::mesh::Selector outputSelector = *bboxMesh.meta_data().get_part("block_1");
    int numNodes = stk::mesh::count_selected_entities(outputSelector & bboxMesh.meta_data().locally_owned_part(), bboxMesh.bulk_data().buckets(stk::topology::NODE_RANK));
    krino::all_reduce_sum(stk::EnvData::parallel_comm(), numNodes);

    output_mesh_and_test_num_nodes(bboxMesh.bulk_data(), outputSelector, numNodes);
  }
}

TEST(OutputUtils, createMeshWithNodeset_outputMeshWithSelectingOnlyNodest_noNodesOutputSinceNotConnectedToAnySelectedElems)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (krino::is_parallel_io_enabled() || stk::parallel_machine_size(comm) == 1)
  {
    krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, comm);

    stk::mesh::Part & nodeset = bboxMesh.meta_data().declare_part_with_topology("MYNODESET", stk::topology::NODE);
    stk::io::put_io_part_attribute(nodeset);

    krino::populate_bounding_box_mesh_and_activate(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 0.3333);

    put_node_with_id_into_nodeset(bboxMesh.bulk_data(), nodeset, 1);

    stk::mesh::Selector outputSelector = nodeset;
    output_mesh_and_test_num_nodes(bboxMesh.bulk_data(), outputSelector, 0);
  }
}

TEST(OutputUtils, decompWithProcBoundarySuchThatOwnedNodeHasNoOwnedElementsSelectedForOutput_noThrowOrHang)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (krino::is_parallel_io_enabled() || stk::parallel_machine_size(comm) == 1)
  {
    krino::BoundingBoxMesh bboxMesh(stk::topology::HEXAHEDRON_8, comm);
    const stk::mesh::Part & block_1 = *bboxMesh.meta_data().get_part("block_1");

    stk::mesh::Part & block_2 = bboxMesh.meta_data().declare_part_with_topology("block_2", stk::topology::HEXAHEDRON_8);
    stk::io::put_io_part_attribute(block_2);

    bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(stk::math::Vector3d{0.,0.,0.}, stk::math::Vector3d{1.,1.,1.}), 0.5);
    bboxMesh.populate_mesh();
    stk::mesh::BulkData & mesh = bboxMesh.bulk_data();

    change_block_for_element_with_id(mesh, block_1, block_2, 8);

    stk::mesh::Selector outputSelector = block_2;
    const krino::ScopedTestFile outputFile("output.e");
    EXPECT_NO_THROW(krino::fix_ownership_and_output_composed_mesh_with_fields(mesh, outputSelector, outputFile.name(), 1, 0.));
  }
}

TEST(OutputUtils, usingMeshFromFileMultipleTimes_makeSureThatDataAddedOntoMetaDataIsRecreatedEachTime)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (krino::is_parallel_io_enabled() || stk::parallel_machine_size(comm) == 1)
  {
    const krino::ScopedTestFile initialMeshFile("mesh.g");
    const std::string & initialMeshName = initialMeshFile.name();
    krino::generate_and_write_bounding_box_mesh(stk::topology::TETRAHEDRON_4, {0.,0.,0.}, {1.,1.,1.}, 0.3333, initialMeshName);

    bool defaultIsTransientFlag = false;
    {
      krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);
      meshFromFile.populate_mesh();

      krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meshFromFile.meta_data());
      defaultIsTransientFlag = cdfemSupport.get_is_transient();
      cdfemSupport.set_is_transient(!defaultIsTransientFlag);

      EXPECT_FALSE(defaultIsTransientFlag == cdfemSupport.get_is_transient());
    }
    {
      krino::MeshFromFile meshFromFile(initialMeshName, MPI_COMM_WORLD);
      meshFromFile.populate_mesh();

      krino::CDFEM_Support & cdfemSupport = krino::CDFEM_Support::get(meshFromFile.meta_data());
      EXPECT_TRUE(defaultIsTransientFlag == cdfemSupport.get_is_transient());
    }
  }
}

namespace
{

// The coordinates of vertex id of the fan written by
// writeFacetMesh_sharedVerticesWrittenOnceWithTheirIds.
stk::math::Vector3d fan_vertex_coordinates(const int64_t id)
{
  return (id == 1) ? stk::math::Vector3d::ZERO : stk::math::Vector3d(1., id - 2., 0.);
}

} // namespace

// A fan of triangles (1, 2+r, 3+r), one owned by each processor r, so vertex 1 is used by every
// processor and vertex 2+r by processors r-1 and r.  The composed file has each vertex once, with
// its id, and the triangles by those ids.
TEST(OutputUtils, writeFacetMesh_sharedVerticesWrittenOnceWithTheirIds)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (!krino::is_parallel_io_enabled() && stk::parallel_machine_size(comm) > 1) return;
  Ioss::Init::Initializer::initialize_ioss();
  const int numProcs = stk::parallel_machine_size(comm);
  const int me = stk::parallel_machine_rank(comm);

  const std::string baseName = "unit_test_facet_mesh";
  const krino::ScopedTestFile outputFile(krino::create_file_name(baseName, 0));
  krino::FacetMesh facets(comm, 3);
  const krino::FacetVertexId first = 2 + me;
  facets.set_owned_facets({{1, fan_vertex_coordinates(1)},
                              {first, fan_vertex_coordinates(first)},
                              {first + 1, fan_vertex_coordinates(first + 1)}},
      {{{1, first, first + 1}}});
  krino::write_facets(facets, baseName, 0);

  if (0 == me)
  {
    Ioss::PropertyManager properties;
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    Ioss::DatabaseIO * db = Ioss::IOFactory::create(
        "exodusII", outputFile.name(), Ioss::READ_MODEL, MPI_COMM_SELF, properties);
    ASSERT_TRUE(db != nullptr && db->ok());
    Ioss::Region io(db, "FacetMeshRegion");

    ASSERT_EQ(numProcs + 2, io.get_property("node_count").get_int());
    Ioss::NodeBlock & nb = *io.get_node_blocks()[0];
    std::vector<int64_t> nodeIds;
    std::vector<double> xyz;
    nb.get_field_data("ids", nodeIds);
    nb.get_field_data("mesh_model_coordinates", xyz);
    std::vector<int64_t> sortedIds(nodeIds);
    std::sort(sortedIds.begin(), sortedIds.end());
    std::vector<int64_t> expectedIds(numProcs + 2);
    std::iota(expectedIds.begin(), expectedIds.end(), 1);
    EXPECT_EQ(expectedIds, sortedIds);
    for (size_t n = 0; n < nodeIds.size(); ++n)
    {
      EXPECT_EQ(fan_vertex_coordinates(nodeIds[n]), stk::math::Vector3d(&xyz[3 * n]))
          << "node " << nodeIds[n];
    }

    Ioss::ElementBlock & eb = *io.get_element_blocks()[0];
    ASSERT_EQ(numProcs, eb.entity_count());
    std::vector<int64_t> elemIds;
    std::vector<int64_t> conn;
    eb.get_field_data("ids", elemIds);
    eb.get_field_data("connectivity", conn);
    for (int e = 0; e < numProcs; ++e)
    {
      const int64_t proc = elemIds[e] - 1;
      EXPECT_EQ((std::vector<int64_t>{1, proc + 2, proc + 3}),
          std::vector<int64_t>(conn.begin() + 3 * e, conn.begin() + 3 * e + 3));
    }
  }
}

// Without shared vertices, processor r's triangle (0,0,r), (1,0,r), (0,1,r) has its own 3 nodes,
// and the nodes and elements are numbered in order of the facets across the processors:
// element r+1 has nodes 3r+1 to 3r+3.
TEST(OutputUtils, writeUnconnectedFacetMesh_eachFacetHasItsOwnNodesNumberedInOrder)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (!krino::is_parallel_io_enabled() && stk::parallel_machine_size(comm) > 1) return;
  Ioss::Init::Initializer::initialize_ioss();
  const int numProcs = stk::parallel_machine_size(comm);
  const int me = stk::parallel_machine_rank(comm);
  auto triangle_vertex = [](const int proc, const int n)
  { return stk::math::Vector3d((n == 1) ? 1. : 0., (n == 2) ? 1. : 0., proc); };

  const std::string baseName = "unit_test_unconnected_facet_mesh";
  const krino::ScopedTestFile outputFile(krino::create_file_name(baseName, 0));
  krino::FacetMesh facets(comm, 3, {.sharedVertices = false});
  facets.set_owned_facets_from_coordinates(
      {triangle_vertex(me, 0), triangle_vertex(me, 1), triangle_vertex(me, 2)});
  krino::write_facets(facets, baseName, 0);

  if (0 == me)
  {
    Ioss::PropertyManager properties;
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    Ioss::DatabaseIO * db = Ioss::IOFactory::create(
        "exodusII", outputFile.name(), Ioss::READ_MODEL, MPI_COMM_SELF, properties);
    ASSERT_TRUE(db != nullptr && db->ok());
    Ioss::Region io(db, "UnconnectedFacetMeshRegion");

    ASSERT_EQ(3 * numProcs, io.get_property("node_count").get_int());
    Ioss::NodeBlock & nb = *io.get_node_blocks()[0];
    std::vector<int64_t> nodeIds;
    std::vector<double> xyz;
    nb.get_field_data("ids", nodeIds);
    nb.get_field_data("mesh_model_coordinates", xyz);
    for (size_t n = 0; n < nodeIds.size(); ++n)
    {
      const int64_t i = nodeIds[n] - 1;
      EXPECT_EQ(triangle_vertex(i / 3, i % 3), stk::math::Vector3d(&xyz[3 * n]))
          << "node " << nodeIds[n];
    }

    Ioss::ElementBlock & eb = *io.get_element_blocks()[0];
    ASSERT_EQ(numProcs, eb.entity_count());
    std::vector<int64_t> elemIds;
    std::vector<int64_t> conn;
    eb.get_field_data("ids", elemIds);
    eb.get_field_data("connectivity", conn);
    for (int e = 0; e < numProcs; ++e)
    {
      const int64_t first = 3 * (elemIds[e] - 1) + 1;
      EXPECT_EQ((std::vector<int64_t>{first, first + 1, first + 2}),
          std::vector<int64_t>(conn.begin() + 3 * e, conn.begin() + 3 * e + 3));
    }
  }
}

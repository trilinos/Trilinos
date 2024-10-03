#include <Akri_AllReduce.hpp>
#include <gtest/gtest.h>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_OutputUtils.hpp>
#include <Ioss_GroupingEntity.h>
#include <Ioss_IOFactory.h>
#include <stk_io/IossBridge.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/EnvData.hpp>

void put_node_with_id_into_nodeset(stk::mesh::BulkData & mesh, stk::mesh::Part & nodeset, const stk::mesh::EntityId id)
{
  const stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, 1);
  mesh.modification_begin();
  if (mesh.is_valid(node) && mesh.bucket(node).owned())
    mesh.change_entity_parts(node, stk::mesh::PartVector{&nodeset});
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
      ThrowRuntimeError("error reading file " << filename);
    }
    std::unique_ptr<Ioss::Region> io = std::make_unique<Ioss::Region>(db, "EXOSurface IC Region");

    numNodes = io->get_property("node_count").get_int();
  }
  krino::all_reduce_sum(stk::EnvData::parallel_comm(), numNodes);
  return numNodes;
}

static void generate_bounding_box_mesh(krino::BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize)
{
  bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), meshSize);
  bboxMesh.set_mesh_structure_type(krino::FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  bboxMesh.populate_mesh();
}

void output_mesh_and_test_num_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const int goldNumOutputNodes)
{
  const std::string filename = "outputTest.e";
  krino::output_composed_mesh_with_fields(mesh, outputSelector, filename, 1, 0.);
  EXPECT_EQ(goldNumOutputNodes, get_num_nodes_in_exodus_file(filename));
}

TEST(OutputUtils, createMesh_outputMeshWithSelectingBlock_allNodesOutput)
{
  krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);
  generate_bounding_box_mesh(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 0.3333);

  stk::mesh::Selector outputSelector = *bboxMesh.meta_data().get_part("block_1");
  int numNodes = stk::mesh::count_selected_entities(outputSelector & bboxMesh.meta_data().locally_owned_part(), bboxMesh.bulk_data().buckets(stk::topology::NODE_RANK));
  krino::all_reduce_sum(stk::EnvData::parallel_comm(), numNodes);

  output_mesh_and_test_num_nodes(bboxMesh.bulk_data(), outputSelector, numNodes);
}

TEST(OutputUtils, createMeshWithNodeset_outputMeshWithSelectingOnlyNodest_noNodesOutputSinceNotConnectedToAnySelectedElems)
{
  krino::BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);

  stk::mesh::Part & nodeset = bboxMesh.meta_data().declare_part_with_topology("MYNODESET", stk::topology::NODE);
  stk::io::put_io_part_attribute(nodeset);

  generate_bounding_box_mesh(bboxMesh, {0.,0.,0.}, {1.,1.,1.}, 0.3333);

  put_node_with_id_into_nodeset(bboxMesh.bulk_data(), nodeset, 1);

  stk::mesh::Selector outputSelector = nodeset;
  output_mesh_and_test_num_nodes(bboxMesh.bulk_data(), outputSelector, 0);
}

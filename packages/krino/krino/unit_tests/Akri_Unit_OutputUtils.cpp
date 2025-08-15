#include <Akri_AllReduce.hpp>
#include <gtest/gtest.h>

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <Ioss_GroupingEntity.h>
#include <Ioss_IOFactory.h>
#include <stk_io/IossBridge.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/environment/EnvData.hpp>

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
      ThrowRuntimeError("error reading file " << filename);
    }
    std::unique_ptr<Ioss::Region> io = std::make_unique<Ioss::Region>(db, "EXOSurface IC Region");

    numNodes = io->get_property("node_count").get_int();
  }
  krino::all_reduce_sum(stk::EnvData::parallel_comm(), numNodes);
  return numNodes;
}

void output_mesh_and_test_num_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const int goldNumOutputNodes)
{
  const std::string filename = "outputTest.e";
  krino::output_composed_mesh_with_fields(mesh, outputSelector, filename, 1, 0.);
  EXPECT_EQ(goldNumOutputNodes, get_num_nodes_in_exodus_file(filename));
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
    EXPECT_NO_THROW(krino::fix_ownership_and_output_composed_mesh_with_fields(mesh, outputSelector, "output.e", 1, 0.));
  }
}

TEST(OutputUtils, usingMeshFromFileMultipleTimes_makeSureThatDataAddedOntoMetaDataIsRecreatedEachTime)
{
  const stk::ParallelMachine comm{MPI_COMM_WORLD};
  if (krino::is_parallel_io_enabled() || stk::parallel_machine_size(comm) == 1)
  {
    const std::string initialMeshName = "mesh.g";
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

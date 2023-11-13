#include "gtest/gtest.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_io/FillMesh.hpp>

namespace
{

TEST(NoUpwardConnectivity, basic_mesh_read_no_aura)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  const unsigned numNodesLocalWithoutAura = numProcs==1 ? 12 : 8;
  const unsigned numElemsLocalWithoutAura = numProcs==1 ? 2 : 1;
  EXPECT_EQ(numNodesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, meta.universal_part()));
  EXPECT_EQ(numElemsLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part()));

  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity node) {
                                   EXPECT_EQ(0u, bulkData.num_elements(node));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::ELEM_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity elem) {
                                   EXPECT_EQ(8u, bulkData.num_nodes(elem));
                                 });
}

TEST(NoUpwardConnectivity, basic_mesh_read_no_aura_create_edges)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  stk::mesh::create_edges(bulk);

  const unsigned numNodesLocalWithoutAura = numProcs==1 ? 12 : 8;
  const unsigned numEdgesLocalWithoutAura = numProcs==1 ? 20 : 12;
  const unsigned numElemsLocalWithoutAura = numProcs==1 ?  2 : 1;
  EXPECT_EQ(numNodesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, meta.universal_part()));
  EXPECT_EQ(numEdgesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::EDGE_RANK, meta.universal_part()));
  EXPECT_EQ(numElemsLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part()));

  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity node) {
                                   EXPECT_EQ(0u, bulkData.num_elements(node));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::EDGE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity edge) {
                                   EXPECT_EQ(0u, bulkData.num_elements(edge));
                                   EXPECT_EQ(2u, bulkData.num_nodes(edge));
                                   const bool onlyDownwardRelations = true;
                                   EXPECT_EQ(bulkData.count_relations(edge,onlyDownwardRelations),
                                             bulkData.count_relations(edge,!onlyDownwardRelations));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::ELEM_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity elem) {
                                   EXPECT_EQ(12u, bulkData.num_edges(elem));
                                   EXPECT_EQ(8u, bulkData.num_nodes(elem));
                                 });
}

TEST(NoUpwardConnectivity, basic_mesh_read_no_aura_create_faces)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  stk::mesh::create_faces(bulk);

  const unsigned numNodesLocalWithoutAura = numProcs==1 ? 12 : 8;
  const unsigned numFacesLocalWithoutAura = numProcs==1 ? 11 : 6;
  const unsigned numElemsLocalWithoutAura = numProcs==1 ?  2 : 1;
  EXPECT_EQ(numNodesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, meta.universal_part()));
  EXPECT_EQ(numFacesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, meta.universal_part()));
  EXPECT_EQ(numElemsLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part()));

  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity node) {
                                   EXPECT_EQ(0u, bulkData.num_elements(node));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::FACE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity face) {
                                   EXPECT_EQ(0u, bulkData.num_elements(face));
                                   EXPECT_EQ(4u, bulkData.num_nodes(face));
                                   const bool onlyDownwardRelations = true;
                                   EXPECT_EQ(bulkData.count_relations(face,onlyDownwardRelations),
                                             bulkData.count_relations(face,!onlyDownwardRelations));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::ELEM_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity elem) {
                                   EXPECT_EQ(6u, bulkData.num_faces(elem));
                                   EXPECT_EQ(8u, bulkData.num_nodes(elem));
                                 });
}

TEST(NoUpwardConnectivity, basic_mesh_read_no_aura_create_edges_and_faces)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  stk::mesh::create_edges(bulk);
  stk::mesh::create_faces(bulk);

  const unsigned numNodesLocalWithoutAura = numProcs==1 ? 12 : 8;
  const unsigned numEdgesLocalWithoutAura = numProcs==1 ? 20 : 12;
  const unsigned numFacesLocalWithoutAura = numProcs==1 ? 11 : 6;
  const unsigned numElemsLocalWithoutAura = numProcs==1 ?  2 : 1;
  EXPECT_EQ(numNodesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, meta.universal_part()));
  EXPECT_EQ(numEdgesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::EDGE_RANK, meta.universal_part()));
  EXPECT_EQ(numFacesLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::FACE_RANK, meta.universal_part()));
  EXPECT_EQ(numElemsLocalWithoutAura, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part()));

  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity node) {
                                   EXPECT_EQ(0u, bulkData.num_elements(node));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::EDGE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity edge) {
                                   EXPECT_EQ(0u, bulkData.num_elements(edge));
                                   EXPECT_EQ(2u, bulkData.num_nodes(edge));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::FACE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity face) {
                                   EXPECT_EQ(0u, bulkData.num_elements(face));
                                   EXPECT_EQ(4u, bulkData.num_nodes(face));
                                 });

  stk::mesh::for_each_entity_run(bulk, stk::topology::ELEM_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity elem) {
                                   EXPECT_EQ(6u, bulkData.num_faces(elem));
                                   EXPECT_EQ(12u, bulkData.num_edges(elem));
                                   EXPECT_EQ(8u, bulkData.num_nodes(elem));
                                 });
}

TEST(NoUpwardConnectivity, elem_graph_not_available)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

  EXPECT_FALSE(bulkPtr->has_face_adjacent_element_graph());
  EXPECT_FALSE(bulkPtr->initialize_face_adjacent_element_graph());
  EXPECT_ANY_THROW(bulkPtr->get_face_adjacent_element_graph());
}

TEST(NoUpwardConnectivity, basic_mesh_read_create_all_sides_throws)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int numProcs = stk::parallel_machine_size(comm);
  if(numProcs > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  EXPECT_ANY_THROW(stk::mesh::create_all_sides(bulk, meta.universal_part()));
}

TEST(NoUpwardConnectivity, basic_mesh_read_with_aura)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) > 2) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(comm);
  builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
  builder.set_upward_connectivity(false);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::io::fill_mesh("generated:1x1x2", bulk);

  const unsigned numNodesLocalWithAura = 12;
  const unsigned numElemsLocalWithAura = 2;
  EXPECT_EQ(numNodesLocalWithAura, stk::mesh::count_entities(bulk, stk::topology::NODE_RANK, meta.universal_part()));
  EXPECT_EQ(numElemsLocalWithAura, stk::mesh::count_entities(bulk, stk::topology::ELEM_RANK, meta.universal_part()));

  stk::mesh::for_each_entity_run(bulk, stk::topology::NODE_RANK, meta.universal_part(),
                                 [](const stk::mesh::BulkData& bulkData, stk::mesh::Entity node) {
                                   EXPECT_EQ(0u, bulkData.num_elements(node));
                                 });
}

}

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "Setup2Block2HexMesh.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_mesh/baseImpl/Visitors.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include <gtest/gtest.h>
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/GetEntities.hpp>       // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <vector>                       // for vector, vector<>::iterator

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::Entity;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::impl::VisitClosureGeneral;
using stk::mesh::impl::VisitClosure;
using stk::mesh::impl::VisitUpwardClosureGeneral;
using stk::mesh::impl::VisitUpwardClosure;
using stk::mesh::impl::VisitAuraClosureGeneral;
using stk::mesh::impl::VisitAuraClosure;
using stk::mesh::impl::VisitUpDownClosure;
using stk::mesh::impl::StoreInVector;
using stk::mesh::impl::StoreInSet;
using stk::mesh::impl::AlwaysVisit;
using stk::mesh::impl::OnlyVisitGhostsOnce;
using stk::mesh::impl::OnlyVisitLocallyOwnedOnce;
using stk::unit_test_util::build_mesh;

TEST ( MeshImplUtils, find_element_edge_ordinal_and_equivalent_nodes )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::BulkData& bulk = *bulkPtr;

  setup2Block2HexMesh(bulk);

  const unsigned numNodesPerEdge = 2;

  //we know that edge 2-6 is edge-ordinal 9 on element 1 and
  //edge-ordinal 8 on element 2
  stk::mesh::EntityId edge_2_6_nodeIds[] = {2, 6};
  stk::mesh::Entity edge_2_6_nodes[numNodesPerEdge];
  edge_2_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[0]);
  edge_2_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[1]);

  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, elem1Id);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, elem2Id);

  stk::mesh::Entity elemEdgeNodes[numNodesPerEdge];
  unsigned elemEdgeOrdinal = 999;

  bool found_it = stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(bulk, elem1, numNodesPerEdge, edge_2_6_nodes,
                                                                                  elemEdgeOrdinal, elemEdgeNodes);

  EXPECT_EQ(true, found_it);
  unsigned expectedElemEdgeOrdinal = 9;
  EXPECT_EQ(expectedElemEdgeOrdinal, elemEdgeOrdinal);
  EXPECT_EQ(edge_2_6_nodes[0], elemEdgeNodes[0]);
  EXPECT_EQ(edge_2_6_nodes[1], elemEdgeNodes[1]);

  found_it = stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(bulk, elem2, numNodesPerEdge, edge_2_6_nodes,
                                                                             elemEdgeOrdinal, elemEdgeNodes);

  EXPECT_EQ(true, found_it);
  expectedElemEdgeOrdinal = 8;
  EXPECT_EQ(expectedElemEdgeOrdinal, elemEdgeOrdinal);
  EXPECT_EQ(edge_2_6_nodes[0], elemEdgeNodes[0]);
  EXPECT_EQ(edge_2_6_nodes[1], elemEdgeNodes[1]);
}

class ClosureFixture
{
public:
  ClosureFixture(MPI_Comm communicator, int num_x, int num_y=1)
  {
    const int spatialDim = 3;
    m_mesh = build_mesh(spatialDim, communicator);
    m_meta = &(m_mesh->mesh_meta_data());
    std::ostringstream oss;
    oss << "generated:" << num_x << "x" << num_y << "x" << m_mesh->parallel_size() << "|sideset:xXyYzZ";
    std::string exodusFileName = stk::unit_test_util::get_option("-i", oss.str());
    stk::io::StkMeshIoBroker exodus_file_reader(communicator);
    exodus_file_reader.set_bulk_data(*m_mesh);
    exodus_file_reader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    exodus_file_reader.create_input_mesh();
    exodus_file_reader.populate_bulk_data();
  }

  ~ClosureFixture()
  {
  }

  void WriteToExodusFile(std::string filename)
  {
    stk::io::StkMeshIoBroker exodus_file_reader(m_mesh->parallel());
    exodus_file_reader.set_bulk_data(*m_mesh);
    int index = exodus_file_reader.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    exodus_file_reader.write_output_mesh(index);
  }

  BulkData & mesh() { return *m_mesh; }
  int prank() { return m_mesh->parallel_rank(); }
  int psize() { return m_mesh->parallel_size(); }

private:
  MetaData * m_meta = nullptr;
  std::shared_ptr<BulkData> m_mesh;
};


std::string PrintEntityVector(EntityVector ev, const BulkData & mesh)
{
  std::ostringstream oss;
  int myRank = mesh.parallel_rank();
  oss << "P" << myRank << " { ";
  for (size_t i=0 ; i<ev.size() ; ++i) {
    EntityKey key = mesh.entity_key(ev[i]);
    oss << key << " ";
  }
  oss << "}\n";
  return oss.str();
}

TEST(MeshImplUtils, visit_closure_trivial)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int spatialDim = 3;
  std::shared_ptr<BulkData> mesh = build_mesh(spatialDim, communicator);

  Entity entity = Entity();
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitClosure(*mesh, entity, siv);
  EXPECT_EQ( 0u, ev.size() );
}

TEST(MeshImplUtils, visit_closure_nominal)

{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();

  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector element_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::ELEMENT_RANK),element_vector);
  ASSERT_TRUE( !element_vector.empty() );
  Entity element = element_vector[0];
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitClosure(mesh, element, siv);
  if (numProcs == 1)
  {
    EXPECT_EQ( 14u, ev.size() );
  }
  else if (myRank == 0 || myRank == numProcs -1)
  {
    EXPECT_EQ( 13u, ev.size() );
  }
  else
  {
    EXPECT_EQ( 12u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_closure_face)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  BulkData & mesh = fix.mesh();

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector face_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::FACE_RANK),face_vector);
  ASSERT_TRUE( !face_vector.empty() );
  Entity face = face_vector[0];
  VisitClosure(mesh, face, siv);
  EXPECT_EQ( 5u, ev.size() );
}

TEST(MeshImplUtils, visit_closure_of_vector)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector element_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::ELEMENT_RANK),element_vector);
  VisitClosure(mesh,element_vector.begin(),element_vector.end(),siv);
  if (numProcs == 1)
  {
    EXPECT_EQ( 24u, ev.size() );
  }
  else if (myRank == 0 || myRank == numProcs-1)
  {
    EXPECT_EQ( 22u, ev.size() );
  }
  else
  {
    EXPECT_EQ( 20u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_closure_of_vector_locally_owned)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector element_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::ELEMENT_RANK),element_vector);
  OnlyVisitLocallyOwnedOnce ovloeo(mesh);
  VisitClosureGeneral(mesh,element_vector.begin(),element_vector.end(),siv,ovloeo);
  if (numProcs == 1)
  {
    EXPECT_EQ( 24u, ev.size() );
  }
  else if (myRank == 0)
  {
    EXPECT_EQ( 22u, ev.size() );
  }
  else if (myRank == numProcs-1)
  {
    EXPECT_EQ( 16u, ev.size() );
  }
  else
  {
    EXPECT_EQ( 14u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_upward_closure_trivial)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int spatialDim = 3;
  std::shared_ptr<BulkData> mesh = build_mesh(spatialDim, communicator);

  Entity entity = Entity();
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitUpwardClosure(*mesh, entity, siv);
  EXPECT_EQ( 0u, ev.size() );
}

TEST(MeshImplUtils, visit_upward_closure_nominal)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 4) { return; }

  Entity node;
  if (myRank == 0) {
    node = mesh.get_entity(stk::topology::NODE_RANK,1);
  }
  else if (myRank == 1) {
    node = mesh.get_entity(stk::topology::NODE_RANK,7);
  }
  else if (myRank == 2) {
    node = mesh.get_entity(stk::topology::NODE_RANK,17);
  }
  else if (myRank == 3) {
    node = mesh.get_entity(stk::topology::NODE_RANK,29);
  }
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitUpwardClosure(mesh, node, siv);
  if (myRank == 0) {
    EXPECT_EQ( 5u, ev.size() );
  }
  else if (myRank == 1) {
    EXPECT_EQ( 7u, ev.size() );
  }
  else if (myRank == 2) {
    EXPECT_EQ( 9u, ev.size() );
  }
  else { // myRank == 3
    EXPECT_EQ( 7u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_upward_closure_face)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  BulkData & mesh = fix.mesh();

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector face_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::FACE_RANK),face_vector);
  ASSERT_TRUE( !face_vector.empty() );
  Entity face = face_vector[0];
  VisitUpwardClosure(mesh, face, siv);
  EXPECT_EQ( 2u, ev.size() );
}

TEST(MeshImplUtils, visit_upward_closure_of_vector)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 4) { return; }

  EntityVector entity_vector;
  if (myRank == 0) {
    entity_vector.push_back(mesh.get_entity(stk::topology::NODE_RANK,1));
    entity_vector.push_back(mesh.get_entity(stk::topology::NODE_RANK,6));
  }
  else if (myRank == 1) {
    entity_vector.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK,3));
  }
  else if (myRank == 2) {
    entity_vector.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK,5));
    entity_vector.push_back(mesh.get_entity(stk::topology::NODE_RANK,15));
  }
  else { // (myRank == 3)
    entity_vector.push_back(mesh.get_entity(stk::topology::NODE_RANK,20));
    entity_vector.push_back(mesh.get_entity(stk::topology::NODE_RANK,23));
  }

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitUpwardClosure(mesh,entity_vector.begin(),entity_vector.end(),siv);
  if (myRank == 0) {
    EXPECT_EQ( 10u, ev.size() );
  }
  else if (myRank == 1) {
    EXPECT_EQ( 1u, ev.size() );
  }
  else if (myRank == 2) {
    EXPECT_EQ( 8u, ev.size() );
  }
  else { // (myRank == 3)
    EXPECT_EQ( 14u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_aura_closure_trivial)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  BulkData & mesh = fix.mesh();

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Entity entity = Entity();
  VisitAuraClosure(mesh,entity,siv);
  EXPECT_EQ( 0u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_element)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs > 4) { return; }

  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector element_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::ELEMENT_RANK),element_vector);
  ASSERT_TRUE( !element_vector.empty() );
  Entity element = element_vector[0];
  VisitAuraClosure(mesh,element,siv);
  if (numProcs == 1)
  {
    EXPECT_EQ( 24u, ev.size() );
  }
  else if (numProcs == 2)
  {
    EXPECT_EQ( 38u, ev.size() );
  }
  else if (myRank == 0 || myRank == numProcs-1)
  {
    EXPECT_EQ( 36u, ev.size() );
  }
  else
  {
    if (numProcs==3) {
      EXPECT_EQ( 52u, ev.size() );
    }
    else {
      EXPECT_EQ( 50u, ev.size() );
    }
  }

}

TEST(MeshImplUtils, visit_aura_closure_of_corner_node)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs > 1) { return; }

  Entity node;
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector node_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::NODE_RANK),node_vector);
  ASSERT_TRUE( !node_vector.empty() );
  node = node_vector[0];
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 14u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_center_node)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs > 1) { return; }

  Entity node;
  Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();
  EntityVector node_vector;
  get_selected_entities(locally_owned_selector, mesh.buckets(stk::topology::NODE_RANK),node_vector);
  ASSERT_TRUE( node_vector.size() > 1 );
  node = node_vector[1];
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 24u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_corner_node_in_2procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs != 2 ) { return; }

  Selector shared_selector = mesh.mesh_meta_data().globally_shared_part();
  EntityVector node_vector;
  get_selected_entities(shared_selector, mesh.buckets(stk::topology::NODE_RANK),node_vector);
  ASSERT_TRUE( !node_vector.empty() );
  Entity node = node_vector[0];
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 22u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_center_node_in_2procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,2);
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs != 2 ) { return; }

  Selector shared_selector = mesh.mesh_meta_data().globally_shared_part();
  EntityVector node_vector;
  get_selected_entities(shared_selector, mesh.buckets(stk::topology::NODE_RANK),node_vector);
  ASSERT_TRUE( node_vector.size() > 1 );
  Entity node = node_vector[1];
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 38u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_side_node_in_3procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4 x 4 x 3
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 3 ) { return; }

  Entity node;
  if (myRank == 0 || myRank == 1) {
    node = mesh.get_entity(stk::topology::NODE_RANK,31);
  }
  else  {
    node = mesh.get_entity(stk::topology::NODE_RANK,56);
  }
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 30u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_center_node_in_3procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4 x 4 x 3
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 3 ) { return; }

  Entity node;
  if (myRank == 0 || myRank == 1) {
    node = mesh.get_entity(stk::topology::NODE_RANK,32);
  }
  else  {
    node = mesh.get_entity(stk::topology::NODE_RANK,57);
  }
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node,siv);
  EXPECT_EQ( 47u, ev.size() );
}

TEST(MeshImplUtils, visit_aura_closure_of_side_element_in_3procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4 x 4 x 3
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 3 ) { return; }

  Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,21);
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,element,siv);
  if (myRank == 0 || myRank == 2) {
    EXPECT_EQ( 64u, ev.size() );
  }
  else {
    EXPECT_EQ( 93u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_aura_closure_of_center_element_in_3procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4 x 4 x 3
  int numProcs = fix.psize();
  const int myRank = fix.prank();
  BulkData & mesh = fix.mesh();
  if (numProcs != 3 ) { return; }

  Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,22);
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,element,siv);
  if (myRank == 0 || myRank == 2) {
    EXPECT_EQ( 87u, ev.size() );
  }
  else  {
    EXPECT_EQ( 127u, ev.size() );
  }
}

TEST(MeshImplUtils, visit_up_down_closure_of_center_element_in_3procs)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4 x 4 x 3
  int numProcs = fix.psize();
  if (numProcs != 3 ) { return; }

  BulkData & mesh = fix.mesh();

  Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,22);
  EntityVector evUpDown;
  StoreInVector<EntityVector> sivUpDown(evUpDown);
  VisitUpDownClosure(mesh,element,sivUpDown);
  EXPECT_EQ( 9u, evUpDown.size() );
}

TEST(MeshImplUtils, visit_aura_closure_vector)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,4,4); // 4x4x1
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs != 1 ) { return; }

  std::vector<Entity> node_vec;
  node_vec.push_back(mesh.get_entity(stk::topology::NODE_RANK,25));
  node_vec.push_back(mesh.get_entity(stk::topology::NODE_RANK,26));
  EntityVector ev;
  StoreInVector<EntityVector> siv(ev);
  VisitAuraClosure(mesh,node_vec.begin(),node_vec.end(),siv);
  EXPECT_EQ( 26u, ev.size() );
}


TEST(MeshImplUtils, visit_aura_closure_vector_ghost)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  ClosureFixture fix(communicator,1,1); // 1x1x4
  const int myRank = fix.prank();
  int numProcs = fix.psize();
  BulkData & mesh = fix.mesh();
  if (numProcs != 4 ) { return; }

  EntityVector ev;
  if (myRank > 0) {
    ev.push_back(mesh.get_entity(stk::topology::NODE_RANK,13));
    ev.push_back(mesh.get_entity(stk::topology::NODE_RANK,14));
    ev.push_back(mesh.get_entity(stk::topology::NODE_RANK,15));
    ev.push_back(mesh.get_entity(stk::topology::NODE_RANK,16));
    ev.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK,3));
    std::set<Entity> es;
    StoreInSet<std::set<Entity> > sis(es);
    OnlyVisitGhostsOnce ovgo(mesh);
    VisitAuraClosureGeneral(mesh,ev.begin(),ev.end(),sis,ovgo);
    if (myRank == 1) {
      EXPECT_EQ( 9u, es.size() );
    }
    else if (myRank == 2) {
      EXPECT_EQ( 19u, es.size() );
    }
    else if (myRank == 3) {
      EXPECT_EQ( 9u, es.size() );
    }
  }
}

TEST(MeshImplUtils, check_for_connected_nodes)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  unsigned spatialDim = 2;
  std::shared_ptr<BulkData> meshPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& meta = meshPtr->mesh_meta_data();
  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::BulkData& mesh = *meshPtr;
  mesh.modification_begin();
  Entity node1, node2, node3, node4;
  if (mesh.parallel_rank() == 0) {
    Entity element = mesh.declare_element(1, stk::mesh::ConstPartVector{&block_1});
    node1 = mesh.declare_node(1);
    node2 = mesh.declare_node(2);
    node3 = mesh.declare_node(3);
    node4 = mesh.declare_node(4);
    //before relations declared, this check should fail
    EXPECT_EQ(-1, stk::mesh::impl::check_for_connected_nodes(mesh));
    mesh.declare_relation(element,node1,0);
    mesh.declare_relation(element,node2,1);
    mesh.declare_relation(element,node3,2);
    mesh.declare_relation(element,node4,3);
  }
  mesh.modification_end();
  //nominal checking here
  EXPECT_EQ(0, stk::mesh::impl::check_for_connected_nodes(mesh));
}

TEST(MeshImplUtils, check_for_connected_nodes_invalid_topology)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  unsigned spatialDim = 2;
  std::shared_ptr<BulkData> meshPtr = build_mesh(spatialDim, communicator);
  meshPtr->modification_begin();
  if (meshPtr->parallel_rank() == 0) {
    meshPtr->declare_element(1);
    //before relations declared, this check should fail
    EXPECT_EQ(-1, stk::mesh::impl::check_for_connected_nodes(*meshPtr));
  }
}

TEST(MeshImplUtils, comm_mesh_very_parallel_consistency_nominal)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::unit_test_util::BulkDataTester mesh(meta, communicator);
  if (mesh.parallel_size() >= 1) {
    mesh.modification_begin();
    Entity node1, node2, node3, node4;
    if (mesh.parallel_rank() == 0) {
      Entity element = mesh.declare_element(1, stk::mesh::ConstPartVector{&block_1});
      node1 = mesh.declare_node(1);
      node2 = mesh.declare_node(2);
      node3 = mesh.declare_node(3);
      node4 = mesh.declare_node(4);
      mesh.declare_relation(element, node1, 0);
      mesh.declare_relation(element, node2, 1);
      mesh.declare_relation(element, node3, 2);
      mesh.declare_relation(element, node4, 3);
    }
    mesh.modification_end();
    std::ostringstream msg ;
    bool is_consistent = mesh.my_comm_mesh_verify_parallel_consistency(msg);
    STK_ThrowErrorMsgIf(!is_consistent, msg.str());
  }
}

TEST(MeshImplUtils, check_no_shared_elements_or_higher_nominal)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  int myRank = stk::parallel_machine_rank(pm);
  if (numProcs >= 2) {
    unsigned spatialDim = 2;
    std::shared_ptr<BulkData> meshPtr = build_mesh(spatialDim, pm);
    stk::mesh::BulkData& mesh = *meshPtr;
    stk::mesh::MetaData& meta = meshPtr->mesh_meta_data();
    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
    mesh.modification_begin();
    Entity node1, node2, node3, node4, node5, node6;
    if (myRank == 0)
    {
      Entity element = mesh.declare_element(1, stk::mesh::ConstPartVector{&block_1});
      node1 = mesh.declare_node(1);
      node2 = mesh.declare_node(2); // shared
      node3 = mesh.declare_node(3); // shared
      node4 = mesh.declare_node(4);
      mesh.declare_relation(element,node1,0);
      mesh.declare_relation(element,node2,1);
      mesh.declare_relation(element,node3,2);
      mesh.declare_relation(element,node4,3);
      mesh.add_node_sharing(node2,1);
      mesh.add_node_sharing(node3,1);
    }
    if (myRank == 1)
    {
      Entity element = mesh.declare_element(2, stk::mesh::ConstPartVector{&block_1});
      node2 = mesh.declare_node(2); // shared
      node5 = mesh.declare_node(5);
      node6 = mesh.declare_node(6);
      node3 = mesh.declare_node(3); // shared
      mesh.declare_relation(element,node2,0);
      mesh.declare_relation(element,node5,1);
      mesh.declare_relation(element,node6,2);
      mesh.declare_relation(element,node3,3);
      mesh.add_node_sharing(node2,0);
      mesh.add_node_sharing(node3,0);
    }
    mesh.modification_end();
    //nominal test of this only, since this is enforced in many places, impossible to create shared elements
    EXPECT_EQ(0, stk::mesh::impl::check_no_shared_elements_or_higher(mesh));
  }
}

void call_get_or_create_face_at_element_side_and_check(stk::mesh::BulkData & mesh, stk::mesh::Entity element, unsigned side_ordinal, unsigned new_face_global_id, stk::mesh::Part & part) {
  mesh.modification_begin();
  stk::mesh::PartVector add_parts(1, &part);
  stk::mesh::Entity new_face = stk::mesh::impl::get_or_create_face_at_element_side(mesh, element, side_ordinal, new_face_global_id, stk::mesh::PartVector(1,&part));
  mesh.modification_end();
  ASSERT_TRUE( mesh.is_valid(new_face) );
  EXPECT_EQ( new_face_global_id, mesh.identifier(new_face) );
  ASSERT_EQ( 2u, mesh.num_elements(new_face));
  EXPECT_TRUE( mesh.bucket(new_face).member(part) );
}

TEST( MeshImplUtils, test_create_face_for_sideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 1) {
    return;
  }
  stk::io::StkMeshIoBroker fixture(pm);

  fixture.add_mesh_database("generated:1x1x2", stk::io::READ_MESH);
  fixture.create_input_mesh();
  stk::topology quad4_topology = stk::topology::QUAD_4;
  stk::mesh::Part & quad4_part = fixture.meta_data().get_topology_root_part(quad4_topology);

  stk::mesh::Part & new_topology_sub_part = fixture.meta_data().declare_part("My Fancy Part",stk::topology::FACE_RANK);
  fixture.meta_data().declare_part_subset(quad4_part, new_topology_sub_part);
  fixture.populate_bulk_data();

  stk::mesh::BulkData & mesh = fixture.bulk_data();

  unsigned elem_global_id = 1;
  stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,elem_global_id);
  unsigned side_ordinal = 5;
  unsigned new_face_global_id = 16;


  stk::mesh::Entity new_face = mesh.get_entity(stk::topology::FACE_RANK, new_face_global_id);
  EXPECT_FALSE( mesh.is_valid(new_face) );
  call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
  call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
  call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,new_topology_sub_part);
}


TEST( MeshImplUtils, test_connect_face_to_other_elements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 1) {
    return;
  }
  stk::io::StkMeshIoBroker fixture(pm);

  fixture.add_mesh_database("generated:1x1x2", stk::io::READ_MESH);
  fixture.create_input_mesh();
  fixture.populate_bulk_data();

  stk::mesh::BulkData & mesh = fixture.bulk_data();

  stk::topology quad4_topology = stk::topology::QUAD_4;
  stk::mesh::Part & quad4_part = fixture.meta_data().get_topology_root_part(quad4_topology);

  unsigned elem_global_id = 1;
  stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,elem_global_id);
  unsigned side_ordinal = 5;
  unsigned new_face_global_id = 16;

  call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
  stk::mesh::Entity new_face = mesh.get_entity(stk::topology::FACE_RANK, new_face_global_id);
  ASSERT_TRUE( mesh.is_valid(new_face) );
  mesh.modification_begin();
  stk::mesh::impl::connect_face_to_other_elements(mesh,new_face,element,side_ordinal);
  mesh.modification_end();

  ASSERT_EQ( 2u, mesh.num_elements(new_face));
  unsigned other_element_global_id = 2;
  stk::mesh::Entity other_element = mesh.get_entity(stk::topology::ELEMENT_RANK, other_element_global_id);
  ASSERT_EQ( 1u, mesh.num_faces(other_element) );
  stk::mesh::Entity attached_face = *mesh.begin_faces(other_element);
  EXPECT_EQ( new_face, attached_face );

  unsigned attached_side_ordinal = *mesh.begin_face_ordinals(other_element);
  const unsigned expected_other_element_side_ordinal = 4;
  EXPECT_EQ( expected_other_element_side_ordinal, attached_side_ordinal );
}


TEST( MeshImplUtils, test_create_shell_status) {
  std::vector<stk::topology> element_topology_vector;
  element_topology_vector.push_back(stk::topology(stk::topology::HEX_8));
  stk::topology original_element_topology = stk::topology::HEX_8;
  std::vector<stk::mesh::impl::ShellStatus> element_shell_status;
  stk::mesh::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
  ASSERT_EQ( 1u, element_shell_status.size() );
  EXPECT_EQ( stk::mesh::impl::NO_SHELLS, element_shell_status[0] );

  element_topology_vector.push_back(stk::topology(stk::topology::HEX_8));
  stk::mesh::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
  ASSERT_EQ( 2u, element_shell_status.size() );
  EXPECT_EQ( stk::mesh::impl::NO_SHELLS, element_shell_status[0] );
  EXPECT_EQ( stk::mesh::impl::NO_SHELLS, element_shell_status[1] );

  element_topology_vector.push_back(stk::topology(stk::topology::SHELL_QUAD_4));
  stk::mesh::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
  ASSERT_EQ( 3u, element_shell_status.size() );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[0] );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[1] );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[2] );

  original_element_topology = stk::topology::SHELL_QUAD_4;
  stk::mesh::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
  ASSERT_EQ( 3u, element_shell_status.size() );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[0] );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[1] );
  EXPECT_EQ( stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[2] );

}

TEST( MeshImplUtils, test_connect_face_to_other_elements_2)
{
  std::vector<int> face_nodes(4);
  face_nodes[0] = 5;
  face_nodes[1] = 6;
  face_nodes[2] = 7;
  face_nodes[3] = 8;

  stk::topology element_side_topology = stk::topology::QUAD_4;

  std::vector<int> element_side_nodes(4);
  element_side_nodes[0] = 5;
  element_side_nodes[1] = 6;
  element_side_nodes[2] = 7;
  element_side_nodes[3] = 8;

  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::NO_SHELLS));
  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

  element_side_nodes[0] = 5;
  element_side_nodes[1] = 8;
  element_side_nodes[2] = 7;
  element_side_nodes[3] = 6;

  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::NO_SHELLS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

  element_side_nodes[0] = 5;
  element_side_nodes[1] = 6;
  element_side_nodes[2] = 7;
  element_side_nodes[3] = 9;

  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::NO_SHELLS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

  face_nodes[0] = 5;
  face_nodes[1] = 8;
  face_nodes[2] = 7;
  face_nodes[3] = 6;

  element_side_nodes[0] = 5;
  element_side_nodes[1] = 6;
  element_side_nodes[2] = 7;
  element_side_nodes[3] = 8;

  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::NO_SHELLS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

  element_side_nodes[0] = 5;
  element_side_nodes[1] = 8;
  element_side_nodes[2] = 7;
  element_side_nodes[3] = 6;

  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::NO_SHELLS));
  EXPECT_TRUE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
  EXPECT_FALSE(stk::mesh::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::mesh::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));
}

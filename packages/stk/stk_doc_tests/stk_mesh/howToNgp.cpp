#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpAtomics.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/stk_config.h>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/getOption.h>

namespace {

using IntDualViewType = Kokkos::DualView<int*, stk::ngp::ExecSpace>;

void set_field_on_device(stk::mesh::BulkData &bulk,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Part &meshPart,
                         stk::mesh::Field<double> &meshField,
                         double fieldVal)
{
  //BEGINNgpSetFieldOnDevice
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  EXPECT_EQ(bulk.mesh_meta_data().spatial_dimension(), ngpMesh.get_spatial_dimension());

  stk::mesh::NgpField<double>& ngpMeshField = stk::mesh::get_updated_ngp_field<double>(meshField);
  EXPECT_EQ(meshField.mesh_meta_data_ordinal(), ngpMeshField.get_ordinal());

  ngpMeshField.clear_sync_state();

  stk::mesh::for_each_entity_run(ngpMesh, rank, meshPart,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   ngpMeshField(entity, 0) = fieldVal;
                                 });

  ngpMeshField.modify_on_device();
  //ENDNgpSetFieldOnDevice
}

template <typename T>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         stk::mesh::Field<T> & stkField,
                         T expectedFieldValue)
{
  //BEGINNgpReadFieldOnHost
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());

  stkField.sync_to_host();

  for(const stk::mesh::Bucket* bptr : buckets) {
    for(stk::mesh::Entity elem : *bptr) {
      const double* fieldData = stk::mesh::field_data(stkField, elem);
      EXPECT_EQ(*fieldData, expectedFieldValue);
    }
  }
  //ENDNgpReadFieldOnHost
}

class NgpHowTo : public stk::unit_test_util::MeshFixture
{
public:
  void setup_test_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    extraPart = &get_meta().declare_part("extraPart");
    stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
    auto &shellQuadField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "myField");
    stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, nullptr);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
        0,2,SHELL_QUAD_4,5,6,7,8";
        stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
  const stk::mesh::Part* extraPart = nullptr;
};

TEST_F(NgpHowTo, loopOverSubsetOfMesh)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
  auto &shellQuadField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "myField");
  stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, nullptr);
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
      0,2,SHELL_QUAD_4,5,6,7,8";
      stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  double fieldVal = 13.0;
  set_field_on_device(get_bulk(), stk::topology::ELEM_RANK, shellQuadPart, shellQuadField, fieldVal);

  shellQuadField.sync_to_host();

  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, shellQuadPart))
  {
    for(stk::mesh::Entity elem : *bucket)
    {
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(shellQuadField, elem));
    }
  }
}

template<typename MeshType>
void test_mesh_up_to_date(stk::mesh::BulkData& bulk, const stk::mesh::Part* extraPart)
{
  //BEGINNgpMeshUpToDate
  MeshType& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  EXPECT_TRUE(ngpMesh.is_up_to_date());

  bulk.modification_begin();
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  bulk.change_entity_parts(node1, stk::mesh::ConstPartVector{extraPart});
  bulk.modification_end();

  EXPECT_FALSE(ngpMesh.is_up_to_date());

  {
    MeshType& newNgpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    EXPECT_TRUE(newNgpMesh.is_up_to_date());
  }

  EXPECT_TRUE(ngpMesh.is_up_to_date());

  //ENDNgpMeshUpToDate
}

TEST_F(NgpHowTo, checkIfUpToDate)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }
  setup_test_mesh();
  test_mesh_up_to_date<stk::mesh::NgpMesh>(get_bulk(), extraPart);
}

template <typename FieldType>
void test_field_on_subset_of_mesh(const stk::mesh::BulkData& bulk, const FieldType& field,
                                  stk::mesh::PartOrdinal partThatHasField,
                                  stk::mesh::PartOrdinal partThatDoesntHaveField)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&field, &bucket, ngpMesh, &partThatHasField, &partThatDoesntHaveField] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           if (bucket.member(partThatHasField)) {
                             STK_NGP_ThrowRequire(field.get_num_components_per_entity(elemIndex) > 0);
                           }
                           if (bucket.member(partThatDoesntHaveField)) {
                             STK_NGP_ThrowRequire(field.get_num_components_per_entity(elemIndex) == 0);
                           }
                         });
                       });
}

TEST_F(NgpHowTo, fieldOnSubsetOfMesh)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
  const stk::mesh::Part &hex8Part = get_meta().get_topology_root_part(stk::topology::HEX_8);

  auto &shellQuadField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "myField");
  stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, nullptr);
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
      0,2,SHELL_QUAD_4,5,6,7,8";
      stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  double fieldVal = 13.0;
  set_field_on_device(get_bulk(), stk::topology::ELEM_RANK, shellQuadPart, shellQuadField, fieldVal);

  shellQuadField.sync_to_host();

  stk::mesh::NgpField<double> & ngpShellFieldAdapter = stk::mesh::get_updated_ngp_field<double>(shellQuadField);

  test_field_on_subset_of_mesh(get_bulk(), ngpShellFieldAdapter,
                               shellQuadPart.mesh_meta_data_ordinal(), hex8Part.mesh_meta_data_ordinal());
}

TEST_F(NgpHowTo, loopOverAllMeshNodes)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  double fieldVal = 13.0;
  set_field_on_device(get_bulk(), stk::topology::NODE_RANK, get_meta().universal_part(), field, fieldVal);

  field.sync_to_host();
  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().universal_part())) {
    for(stk::mesh::Entity node : *bucket) {
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(field, node));
    }
  }
}

TEST_F(NgpHowTo, loopOverMeshFaces)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &facePart = get_meta().declare_part("facePart", stk::topology::FACE_RANK);
  auto &field = get_meta().declare_field<double>(stk::topology::FACE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, facePart, nullptr);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {&facePart});

  double fieldVal = 13.0;
  set_field_on_device(get_bulk(), stk::topology::FACE_RANK, facePart, field, fieldVal);

  field.sync_to_host();

  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::FACE_RANK, get_meta().universal_part())) {
    for(stk::mesh::Entity node : *bucket) {
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(field, node));
    }
  }
}

void run_connected_node_test(const stk::mesh::BulkData& bulk)
{
  //BEGINNgpMeshConnectivity
  stk::mesh::Entity elem1_host = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::ConnectedEntities elem1_nodes_host = bulk.get_connected_entities(elem1_host, stk::topology::NODE_RANK);
  stk::mesh::Entity node0_host = elem1_nodes_host[0];
  stk::mesh::Entity node7_host = elem1_nodes_host[7];

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::mesh::Selector allElems = bulk.mesh_meta_data().universal_part();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, allElems,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elemIndex) {
      stk::mesh::Entity elem = ngpMesh.get_entity(stk::topology::ELEM_RANK, elemIndex);
      stk::mesh::EntityId elemId = ngpMesh.identifier(elem);
      if (elemId == 1) {
        STK_NGP_ThrowRequire(elem == elem1_host);
        const stk::mesh::NgpMesh::BucketType& ngpBucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elemIndex.bucket_id);
        STK_NGP_ThrowRequire(ngpBucket.topology() == stk::topology::HEX_8);

        stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);
        STK_NGP_ThrowRequire(nodes.size() == ngpBucket.topology().num_nodes());
        STK_NGP_ThrowRequire(node0_host == nodes[0]);
        STK_NGP_ThrowRequire(node7_host == nodes[7]);

        stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[0]);
        stk::mesh::NgpMesh::ConnectedEntities node0_elems = ngpMesh.get_elements(stk::topology::NODE_RANK, nodeIndex);
        STK_NGP_ThrowRequire(1 == node0_elems.size());
        STK_NGP_ThrowRequire(node0_elems[0] == elem);
      }
    });
  //ENDNgpMeshConnectivity
}

TEST_F(NgpHowTo, loopOverElemNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  run_connected_node_test(get_bulk());
}

TEST_F(NgpHowTo, loopOverElemNodes_bucketCapacity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  const unsigned bucketCapacity = 8;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity, bucketCapacity);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  run_connected_node_test(get_bulk());
}

void run_id_test(const stk::mesh::BulkData& bulk)
{
  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();
                         NGP_EXPECT_EQ(1u, numElems);

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           NGP_EXPECT_EQ(1u, ngpMesh.identifier(elem));
                           NGP_EXPECT_EQ(1u, ngpMesh.entity_key(elem).id());
                           NGP_EXPECT_EQ(stk::topology::ELEM_RANK, ngpMesh.entity_rank(elem));
                           NGP_EXPECT_EQ(stk::topology::ELEM_RANK, ngpMesh.entity_key(elem).rank());

                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);
                           for (unsigned inode = 0; inode < nodes.size(); ++inode) {
                             NGP_EXPECT_EQ(inode+1, ngpMesh.identifier(nodes[inode]));
                             NGP_EXPECT_EQ(inode+1, ngpMesh.entity_key(nodes[inode]).id());
                             NGP_EXPECT_EQ(stk::topology::NODE_RANK, ngpMesh.entity_rank(nodes[inode]));
                             NGP_EXPECT_EQ(stk::topology::NODE_RANK, ngpMesh.entity_key(nodes[inode]).rank());
                           }
                         });
                       });
}

NGP_TEST_F(NgpHowTo, checkElemNodeIds)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  run_id_test(get_bulk());
}

void run_connected_face_test(const stk::mesh::BulkData& bulk)
{
  stk::topology elemTopo = stk::topology::HEX_8;

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
  EXPECT_EQ(1u, elems.size());
  stk::mesh::Entity face0 = bulk.begin_faces(elems[0])[0];
  stk::mesh::Entity face1 = bulk.begin_faces(elems[0])[1];
  stk::mesh::Entity face2 = bulk.begin_faces(elems[0])[2];
  stk::mesh::Entity face3 = bulk.begin_faces(elems[0])[3];
  stk::mesh::Entity face4 = bulk.begin_faces(elems[0])[4];
  stk::mesh::Entity face5 = bulk.begin_faces(elems[0])[5];

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::ConnectedEntities faces = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex);
                           stk::topology bucketTopo = bucket.topology();
                           STK_NGP_ThrowRequire(elemTopo == bucketTopo);
                           STK_NGP_ThrowRequire(faces.size() == bucketTopo.num_faces());
                           STK_NGP_ThrowRequire(face0 == faces[0]);
                           STK_NGP_ThrowRequire(face1 == faces[1]);
                           STK_NGP_ThrowRequire(face2 == faces[2]);
                           STK_NGP_ThrowRequire(face3 == faces[3]);
                           STK_NGP_ThrowRequire(face4 == faces[4]);
                           STK_NGP_ThrowRequire(face5 == faces[5]);

                           stk::mesh::NgpMesh::ConnectedEntities edges = ngpMesh.get_edges(stk::topology::ELEM_RANK, elemIndex);
                           STK_NGP_ThrowRequire(0 == edges.size());

                           stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
                           stk::mesh::NgpMesh::ConnectedEntities face0_elems = ngpMesh.get_elements(stk::topology::FACE_RANK,
                           faceIndex);
                           STK_NGP_ThrowRequire(1 == face0_elems.size());
                           STK_NGP_ThrowRequire(face0_elems[0] == elem);
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFaces)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    GTEST_SKIP();
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1|sideset:xXyYzZ", get_bulk());

  run_connected_face_test(get_bulk());
}

void add_constraint_on_nodes_1_thru_4(stk::mesh::BulkData& bulk,
                                      stk::mesh::EntityId constraintId)
{
  bulk.modification_begin();
  stk::mesh::Entity constraintEntity = bulk.declare_constraint(constraintId);
  for(stk::mesh::EntityId nodeId = 1; nodeId <= 4; ++nodeId) {
    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
    stk::mesh::ConnectivityOrdinal ord = static_cast<stk::mesh::ConnectivityOrdinal>(nodeId-1);
    bulk.declare_relation(constraintEntity, node, ord);
  }
  bulk.modification_end();
}

void run_constraint_node_test(const stk::mesh::BulkData& bulk,
                              stk::mesh::EntityId constraintId)
{
  stk::mesh::Entity constraint = bulk.get_entity(stk::topology::CONSTRAINT_RANK, constraintId);
  ASSERT_TRUE(bulk.is_valid(constraint));
  stk::mesh::Entity node0 = bulk.begin_nodes(constraint)[0];
  stk::mesh::Entity node1 = bulk.begin_nodes(constraint)[1];
  stk::mesh::Entity node2 = bulk.begin_nodes(constraint)[2];
  stk::mesh::Entity node3 = bulk.begin_nodes(constraint)[3];
  ASSERT_TRUE(bulk.is_valid(node0));
  ASSERT_TRUE(bulk.is_valid(node1));
  ASSERT_TRUE(bulk.is_valid(node2));
  ASSERT_TRUE(bulk.is_valid(node3));

  const unsigned expectedNumNodes = bulk.num_nodes(constraint);

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
                       KOKKOS_LAMBDA(const unsigned& i) {
                         stk::mesh::FastMeshIndex constraintMeshIndex = ngpMesh.fast_mesh_index(constraint);

                         stk::mesh::NgpMesh::ConnectedEntities ngpNodes = ngpMesh.get_connected_entities(stk::topology::CONSTRAINT_RANK, constraintMeshIndex, stk::topology::NODE_RANK);
                         NGP_EXPECT_EQ(expectedNumNodes, ngpNodes.size());
                         NGP_EXPECT_EQ(node0, ngpNodes[0]);
                         NGP_EXPECT_EQ(node1, ngpNodes[1]);
                         NGP_EXPECT_EQ(node2, ngpNodes[2]);
                         NGP_EXPECT_EQ(node3, ngpNodes[3]);

                         stk::mesh::NgpMesh::ConnectedOrdinals ngpOrdinals = ngpMesh.get_connected_ordinals(stk::topology::CONSTRAINT_RANK, constraintMeshIndex, stk::topology::NODE_RANK);
                         NGP_EXPECT_EQ(4u, ngpOrdinals.size());
                         NGP_EXPECT_EQ(0, ngpOrdinals[0]);
                         NGP_EXPECT_EQ(1, ngpOrdinals[1]);
                         NGP_EXPECT_EQ(2, ngpOrdinals[2]);
                         NGP_EXPECT_EQ(3, ngpOrdinals[3]);
                       }
                       );
}

class NgpHowToConstraint : public stk::unit_test_util::MeshFixture
{
public:
  NgpHowToConstraint() : MeshFixture(3, {"node", "edge", "face", "elem", "constraint"})
  {}
};

TEST_F(NgpHowToConstraint, checkNodalConnectivity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    GTEST_SKIP();
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());
  stk::mesh::EntityId constraintId = 1;
  add_constraint_on_nodes_1_thru_4(get_bulk(), constraintId);

  run_constraint_node_test(get_bulk(), constraintId);
}

void run_connected_face_ordinal_test(const stk::mesh::BulkData& bulk)
{
  stk::topology elemTopo = stk::topology::HEX_8;

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
  EXPECT_EQ(1u, elems.size());
  stk::mesh::ConnectivityOrdinal ordinal0 = bulk.begin_face_ordinals(elems[0])[0];
  stk::mesh::ConnectivityOrdinal ordinal1 = bulk.begin_face_ordinals(elems[0])[1];
  stk::mesh::ConnectivityOrdinal ordinal2 = bulk.begin_face_ordinals(elems[0])[2];
  stk::mesh::ConnectivityOrdinal ordinal3 = bulk.begin_face_ordinals(elems[0])[3];
  stk::mesh::ConnectivityOrdinal ordinal4 = bulk.begin_face_ordinals(elems[0])[4];
  stk::mesh::ConnectivityOrdinal ordinal5 = bulk.begin_face_ordinals(elems[0])[5];

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::ConnectedOrdinals ordinals = ngpMesh.get_face_ordinals(stk::topology::ELEM_RANK,
                           elemIndex);
                           stk::topology bucketTopo = bucket.topology();
                           STK_NGP_ThrowRequire(elemTopo == bucketTopo);
                           STK_NGP_ThrowRequire(ordinals.size() == bucketTopo.num_faces());
                           STK_NGP_ThrowRequire(ordinal0 == ordinals[0]);
                           STK_NGP_ThrowRequire(ordinal1 == ordinals[1]);
                           STK_NGP_ThrowRequire(ordinal2 == ordinals[2]);
                           STK_NGP_ThrowRequire(ordinal3 == ordinals[3]);
                           STK_NGP_ThrowRequire(ordinal4 == ordinals[4]);
                           STK_NGP_ThrowRequire(ordinal5 == ordinals[5]);

                           stk::mesh::NgpMesh::ConnectedOrdinals edgeOrdinals = ngpMesh.get_edge_ordinals(stk::topology::ELEM_RANK,
                           elemIndex);
                           STK_NGP_ThrowRequire(0 == edgeOrdinals.size());

                           stk::mesh::Entity face = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex)[0];
                           stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(face);
                           stk::mesh::NgpMesh::ConnectedOrdinals face0_elemOrdinals = ngpMesh.get_element_ordinals(stk::topology::FACE_RANK, faceIndex);
                           STK_NGP_ThrowRequire(1 == face0_elemOrdinals.size());
                           STK_NGP_ThrowRequire(face0_elemOrdinals[0] == ordinal0);
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFaceOrdinals)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1|sideset:xXyYzZ", get_bulk());

  run_connected_face_ordinal_test(get_bulk());
}

constexpr unsigned numFacesPerHex = 6;
constexpr unsigned numEdgesPerHex = 12;
constexpr unsigned numHexes = 2;

struct PermutationData {
  stk::mesh::Permutation elemFacePermutations[numHexes][numFacesPerHex];
  stk::mesh::Permutation elemEdgePermutations[numHexes][numEdgesPerHex];
};

void run_connected_face_permutation_test(const stk::mesh::BulkData& bulk)
{
  stk::topology elemTopo = stk::topology::HEX_8;

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
  EXPECT_EQ(2u, elems.size());
  PermutationData p;
  for(unsigned i=0; i<numHexes; ++i) {
    const stk::mesh::Permutation* perms = bulk.begin_face_permutations(elems[i]);
    std::copy(perms, perms+numFacesPerHex, p.elemFacePermutations[i]);
    perms = bulk.begin_edge_permutations(elems[i]);
    std::copy(perms, perms+numEdgesPerHex, p.elemEdgePermutations[i]);
  }

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::Permutations facePermutations = ngpMesh.get_face_permutations(stk::topology::ELEM_RANK, elemIndex);
                           stk::mesh::NgpMesh::Permutations edgePermutations = ngpMesh.get_edge_permutations(stk::topology::ELEM_RANK, elemIndex);
                           stk::topology bucketTopo = bucket.topology();
                           STK_NGP_ThrowRequire(elemTopo == bucketTopo);
                           STK_NGP_ThrowRequire(facePermutations.size() == bucketTopo.num_faces());
                           STK_NGP_ThrowRequire(edgePermutations.size() == bucketTopo.num_edges());

                           for(unsigned j=0; j<numFacesPerHex; ++j) {
                             STK_NGP_ThrowRequire(p.elemFacePermutations[i][j] == facePermutations[j]);
                           }
                           for(unsigned j=0; j<numEdgesPerHex; ++j) {
                             STK_NGP_ThrowRequire(p.elemEdgePermutations[i][j] == edgePermutations[j]);
                           }
                           if (i == 0) {
                             stk::mesh::Entity face = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex)[0];
                             stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(face);
                             stk::mesh::NgpMesh::Permutations face0_elemPermutations = ngpMesh.get_element_permutations(stk::topology::FACE_RANK, faceIndex);
                             STK_NGP_ThrowRequire(1 == face0_elemPermutations.size());
                             STK_NGP_ThrowRequire(face0_elemPermutations[0] == p.elemFacePermutations[0][0]);
                           }
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFacePermutations)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part& sidePart = get_meta().declare_part("SidePart");
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x2|sideset:xXyYzZ", get_bulk());
  stk::mesh::create_edges(get_bulk());
  bool connectFacesToEdges = true;
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {&sidePart}, connectFacesToEdges);

  run_connected_face_permutation_test(get_bulk());
}

void run_another_connected_face_test(const stk::mesh::BulkData& bulk)
{
  stk::topology elemTopo = stk::topology::HEX_8;

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
  EXPECT_EQ(4u, elems.size());

  unsigned numResults = 2;
  auto resultDevice = Kokkos::View<int*, stk::ngp::ExecSpace>("device", numResults);
  auto resultHost = Kokkos::create_mirror_view(resultDevice);
  Kokkos::deep_copy(resultHost, 0);
  Kokkos::deep_copy(resultDevice, 0);

  enum {ELEM_FACE_CHECK = 0, FACE_NODE_CHECK = 1};

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);
  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();
                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&resultDevice, &bucket, ngpMesh, &elemTopo] (const int& i)
                         {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::ConnectedEntities faces = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex);
                           stk::topology bucketTopo = bucket.topology();
                           STK_NGP_ThrowRequire(elemTopo == bucketTopo);
                           int numFaces = faces.size();
                           stk::mesh::atomic_add(&resultDevice(ELEM_FACE_CHECK), numFaces);

                           if (numFaces == 1)
                           {
                             stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
                             stk::mesh::NgpMesh::ConnectedEntities faceNodes = ngpMesh.get_connected_entities(stk::topology::FACE_RANK, faceIndex, stk::topology::NODE_RANK);

                             unsigned faceOrdinal = 0; //when we add ordinals we can fix this. But for a hex all faces have the same topology anyway...
                             stk::topology faceTopo = elemTopo.face_topology(faceOrdinal);

                             STK_NGP_ThrowRequire(faceNodes.size() == faceTopo.num_nodes());
                             stk::mesh::atomic_add(&resultDevice(FACE_NODE_CHECK), 1);
                           }

                         });
                       });

  Kokkos::deep_copy(resultHost, resultDevice);
  EXPECT_EQ(2, resultHost(ELEM_FACE_CHECK)); //expected 2 elements that had faces
  EXPECT_EQ(2, resultHost(FACE_NODE_CHECK)); //expected 2 faces that had nodes
}

TEST_F(NgpHowTo, anotherElemFacesTest)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x4|sideset:zZ", get_bulk());
  run_another_connected_face_test(get_bulk());
}

void test_ngp_mesh_construction(const stk::mesh::BulkData& bulk)
{
  size_t numHostElemBuckets = bulk.buckets(stk::topology::ELEM_RANK).size();
  std::vector<size_t> counts(bulk.mesh_meta_data().entity_rank_count(), 0);
  stk::mesh::count_entities(bulk.mesh_meta_data().locally_owned_part(), bulk, counts);
  size_t numElements = counts[stk::topology::ELEM_RANK];

  double startTime = stk::wall_time();

  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  double elapsedTime = stk::wall_time() - startTime;
  std::cout << "Time to construct stk::mesh::NgpMesh with "<<numElements<<" elements: "<<elapsedTime << std::endl;

  size_t numDeviceElemBuckets = ngpMesh.num_buckets(stk::topology::ELEM_RANK);

  EXPECT_EQ(numHostElemBuckets, numDeviceElemBuckets);
}

TEST_F(NgpHowTo, ngpMeshConstruction)
{
  std::string exodusFileName = stk::unit_test_util::get_option("-mesh", "generated:20x20x20|sideset:xXyYzZ");

  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1 && exodusFileName == "generated:20x20x20|sideset:xXyYzZ") {
    std::cout<<"NgpHowTo.ngpMeshConstruction Only runs in parallel if user specified a mesh." << std::endl;
    return;
  }

  std::cout << "Using mesh: " << exodusFileName << std::endl;

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::MetaData& meta = get_meta();
  stk::mesh::Part& boundaryPart = meta.declare_part("boundary");

  stk::io::fill_mesh(exodusFileName, get_bulk());

  stk::mesh::create_all_block_boundary_sides(get_bulk(), meta.universal_part(), {&boundaryPart});

  test_ngp_mesh_construction(get_bulk());
}

unsigned count_num_elems(stk::mesh::NgpMesh ngpMesh,
                         stk::mesh::NgpField<int> ngpField,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Part &part)
{
  Kokkos::View<unsigned *, stk::ngp::MemSpace> numElems("numElems", 1);
  stk::mesh::for_each_entity_run(ngpMesh, rank, part,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   unsigned fieldValue = static_cast<unsigned>(ngpField(entity, 0));
                                   Kokkos::atomic_add(&numElems(0), fieldValue);
                                 });
  Kokkos::View<unsigned *, stk::ngp::MemSpace>::HostMirror numElemsHost =
      Kokkos::create_mirror_view(numElems);
  Kokkos::deep_copy(numElemsHost, numElems);
  return numElemsHost(0);
}

void set_num_elems_in_field_on_device_and_copy_back(stk::mesh::BulkData &bulk,
                                                    stk::mesh::Part &part,
                                                    stk::mesh::Field<int> &field)
{
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  unsigned numElems = count_num_elems(ngpMesh, ngpField, field.entity_rank(), part);
  stk::mesh::for_each_entity_run(ngpMesh, field.entity_rank(), part,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   ngpField(entity, 0) = numElems;
                                 });
  ngpField.modify_on_device();
  ngpField.sync_to_host();
}

TEST_F(NgpHowTo, exerciseAura)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  auto &field = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "myField");
  int init = 1;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  std::string meshDesc;
  if (get_parallel_size() == 1) {
    meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
               "0,2,HEX_8,5,6,7,8,9,10,11,12";
  }
  else {
    meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
               "1,2,HEX_8,5,6,7,8,9,10,11,12";
  }
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  set_num_elems_in_field_on_device_and_copy_back(get_bulk(), get_meta().universal_part(), field);

  int expectedNumElemsPerProc = 2;
  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().universal_part()))
    for(stk::mesh::Entity elem : *bucket)
      EXPECT_EQ(expectedNumElemsPerProc, *stk::mesh::field_data(field, elem));
}

template <typename DataType>
stk::mesh::Field<DataType> &create_field_with_num_states_and_init(stk::mesh::MetaData &meta, const std::string & fieldName, int numStates, DataType init)
{
  auto &field = meta.declare_field<DataType>(stk::topology::ELEM_RANK, fieldName, numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &init);
  return field;
}

template <typename DataType>
stk::mesh::Field<DataType> &create_vector_field_with_num_states_and_init(stk::mesh::MetaData &meta,
                                                                         const std::string & fieldName,
                                                                         int numStates,
                                                                         int fieldDimension,
                                                                         DataType* init)
{
  auto &field = meta.declare_field<DataType>(stk::topology::ELEM_RANK, fieldName, numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), fieldDimension, init);
  return field;
}

stk::mesh::Field<int> &create_field(stk::mesh::MetaData &meta)
{
  int numStates = 1;
  int initialValue = -1;
  return create_field_with_num_states_and_init(meta, "myField", numStates, initialValue);
}

TEST_F(NgpHowTo, setAllScalarFieldValues)
{
  int numStates = 2;
  double initialValue = 0.0;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::Field<double> &stkField = create_field_with_num_states_and_init<double>(get_meta(), "myField", numStates, initialValue);
  stk::io::fill_mesh("generated:1x1x4", get_bulk());
  stkField.sync_to_device();

  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

  double fieldVal = 1.0;
  ngpField.set_all(ngpMesh, fieldVal);

  double expectedSum = stk::mesh::count_selected_entities(get_meta().locally_owned_part(),
                                                          get_bulk().buckets(stk::topology::ELEM_RANK));
  double sum = stk::mesh::get_field_sum(ngpMesh, ngpField, get_meta().locally_owned_part());
  EXPECT_NEAR(expectedSum, sum, 1e-14);
}

TEST_F(NgpHowTo, setAllVectorFieldValues)
{
  int numStates = 2;
  double initialValue[] = {0.0, 0.0, 0.0};
  int fieldDimension = 3;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  auto& stkField = create_vector_field_with_num_states_and_init<double>(get_meta(),
                                                                        "myField",
                                                                        numStates,
                                                                        fieldDimension,
                                                                        initialValue);
  stk::io::fill_mesh("generated:1x1x4", get_bulk());
  stkField.sync_to_device();

  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());

  double fieldVal = 1.0;
  ngpField.set_all(ngpMesh, fieldVal);

  double expectedSum = fieldDimension*stk::mesh::count_selected_entities(get_meta().locally_owned_part(),
                                                                         get_bulk().buckets(stk::topology::ELEM_RANK));
  double sum = stk::mesh::get_field_sum(ngpMesh, ngpField, get_meta().locally_owned_part());
  EXPECT_NEAR(expectedSum, sum, 1e-14);
}

struct CheckInitialFieldSizeAndValues {
  CheckInitialFieldSizeAndValues(const stk::mesh::NgpField<double>& _ngpField, double _tolerance)
    : ngpField(_ngpField), tolerance(_tolerance)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& meshIdx) const
  {
    NGP_EXPECT_NEAR(1.0, ngpField(meshIdx, 0), tolerance);
    NGP_EXPECT_NEAR(2.0, ngpField(meshIdx, 1), tolerance);
    NGP_EXPECT_NEAR(3.0, ngpField(meshIdx, 2), tolerance);

    stk::mesh::EntityFieldData<double> vals = ngpField(meshIdx);
    NGP_EXPECT_EQ(vals.size(), ngpField.get_num_components_per_entity(meshIdx));
    NGP_EXPECT_NEAR(1.0, vals[0], tolerance);
    NGP_EXPECT_NEAR(2.0, vals[1], tolerance);
    NGP_EXPECT_NEAR(3.0, vals[2], tolerance);
    vals[0] += 2.0;
    vals[1] += 3.0;
    vals[2] += 4.0;
  }

private:
  stk::mesh::NgpField<double> ngpField;
  double tolerance;
};

struct CheckUpdatedFieldSizeAndValues {
  CheckUpdatedFieldSizeAndValues(const stk::mesh::NgpField<double>& _ngpField, double _tolerance)
    : ngpField(_ngpField), tolerance(_tolerance)
  {
  }

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& meshIdx) const
  {
    const stk::mesh::EntityFieldData<double> vals = ngpField(meshIdx);

    NGP_EXPECT_NEAR(3.0, vals[0], tolerance);
    NGP_EXPECT_NEAR(5.0, vals[1], tolerance);
    NGP_EXPECT_NEAR(7.0, vals[2], tolerance);
  }

private:
  stk::mesh::NgpField<double> ngpField;
  double tolerance;
};

void test_vector_field_size_and_values(stk::mesh::BulkData& bulk,
                                       const stk::mesh::FieldBase& stkField)
{
  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  const double tol = 1.e-12;
  const stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part("block_1");
  CheckInitialFieldSizeAndValues checkInitialFieldSizeAndValues(ngpField, tol);
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), block1, checkInitialFieldSizeAndValues);

  CheckUpdatedFieldSizeAndValues checkUpdatedFieldSizeAndValues(ngpField, tol);
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), block1, checkUpdatedFieldSizeAndValues);
}
NGP_TEST_F(NgpHowTo, accessVectorFieldValues)
{
  int numStates = 1;
  double initialValue[] = {1.0, 2.0, 3.0};
  int fieldDimension = 3;
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  auto& stkField = create_vector_field_with_num_states_and_init<double>(get_meta(),
                                                                        "myField",
                                                                        numStates,
                                                                        fieldDimension,
                                                                        initialValue);
  stk::io::fill_mesh("generated:1x1x4", get_bulk());

  test_vector_field_size_and_values(get_bulk(), stkField);
}


class NgpReduceHowTo : public stk::unit_test_util::MeshFixture
{
protected:
  NgpReduceHowTo()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    elemField = &create_field(get_meta());
    stk::io::fill_mesh("generated:1x2x4", get_bulk());
    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);
    for(stk::mesh::Entity elem : elems)
    {
      int *fieldData = stk::mesh::field_data(*elemField, elem);
      fieldData[0] = get_bulk().identifier(elem);
    }
  }
  int get_num_elems()
  {
    return stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
  }
  stk::mesh::Field<int> *elemField;
};

TEST_F(NgpReduceHowTo, getMinFieldValue)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  const int lowestElemID = 1;
  int expectedMinVal = lowestElemID;
  int actualMinVal = stk::mesh::get_field_min(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedMinVal, actualMinVal);
}

TEST_F(NgpReduceHowTo, getMaxFieldValue)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  const int highestElemID = 8;
  int expectedMaxVal = highestElemID;
  int actualMaxVal = stk::mesh::get_field_max(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedMaxVal, actualMaxVal);
}

TEST_F(NgpReduceHowTo, getMaxFieldValueComponent)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  auto * coordField = get_meta().coordinate_field();
  auto & ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*coordField);

  const double expectedMaxCoord = 4.0;
  const double expectedMaxXCoord = 1.0;
  const double expectedMaxYCoord = 2.0;
  double maxVal;
  double maxXVal;
  double maxYVal;
  Kokkos::Max<double> maxValReduction(maxVal);
  Kokkos::Max<double> maxValXReduction(maxXVal);
  Kokkos::Max<double> maxValYReduction(maxYVal);
  stk::mesh::get_field_reduction
      (ngpMesh, ngpCoordField, get_bulk().mesh_meta_data().universal_part(), maxValReduction);
  const int xComponent = 0;
  const int yComponent = 1;
  stk::mesh::get_field_reduction
      (ngpMesh, ngpCoordField, get_bulk().mesh_meta_data().universal_part(), maxValXReduction, xComponent);
  stk::mesh::get_field_reduction
      (ngpMesh, ngpCoordField, get_bulk().mesh_meta_data().universal_part(), maxValYReduction, yComponent);

  EXPECT_EQ(expectedMaxCoord, maxVal);
  EXPECT_EQ(expectedMaxXCoord, maxXVal);
  EXPECT_EQ(expectedMaxYCoord, maxYVal);
}

TEST_F(NgpReduceHowTo, getSumFieldValue)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int numElems = get_num_elems();
  int expectedSum = numElems*(numElems+1)/2;
  int sum_val = stk::mesh::get_field_sum(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedSum, sum_val);
}
TEST_F(NgpReduceHowTo, minMaxPairWiseReduction)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  Kokkos::MinMaxScalar<int> minMaxVal;
  Kokkos::MinMax<int> minMax(minMaxVal);
  stk::mesh::get_field_reduction
      (ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part(), minMax);
  EXPECT_EQ(1, minMaxVal.min_val);
  EXPECT_EQ(get_num_elems(), minMaxVal.max_val);
}
TEST_F(NgpReduceHowTo, minLocReduction)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int expectedMin = 1;
  stk::mesh::EntityId expectedMinLoc = 1;
  Kokkos::ValLocScalar<int,stk::mesh::EntityId> minLocVal;
  Kokkos::MinLoc<int,stk::mesh::EntityId> minLoc (minLocVal);
  stk::mesh::get_field_reduction
      (ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part(), minLoc);
  EXPECT_EQ(expectedMin, minLocVal.val);
  EXPECT_EQ(expectedMinLoc, minLocVal.loc);
}
TEST_F(NgpReduceHowTo, minMaxLocReduction)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(get_bulk());
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int expectedMin = 1;
  stk::mesh::EntityId expectedMinLoc = 1;
  int expectedMax = get_num_elems();
  stk::mesh::EntityId expectedMaxLoc = 8;
  Kokkos::MinMaxLocScalar<int,stk::mesh::EntityId> minMaxLocVal;
  Kokkos::MinMaxLoc<int,stk::mesh::EntityId> minMaxLoc (minMaxLocVal);
  stk::mesh::get_field_reduction (ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part(), minMaxLoc);
  EXPECT_EQ(expectedMin, minMaxLocVal.min_val);
  EXPECT_EQ(expectedMinLoc, minMaxLocVal.min_loc);
  EXPECT_EQ(expectedMax, minMaxLocVal.max_val);
  EXPECT_EQ(expectedMaxLoc, minMaxLocVal.max_loc);
}

template <typename T>
void fill_field_on_device(stk::mesh::BulkData & bulk,
                          stk::mesh::Field<T> & stkField,
                          T fieldValue)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK,
                                 bulk.mesh_meta_data().locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   ngpField(entity, 0) = fieldValue;
                                 });
  ngpField.modify_on_device();
}

template <typename T>
void check_field_on_device(stk::mesh::BulkData & bulk,
                           stk::mesh::Field<T> & stkField,
                           T expectedFieldValue)
{
  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK,
                                 bulk.mesh_meta_data().locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   NGP_EXPECT_EQ(ngpField(entity, 0), expectedFieldValue);
                                 });
}

NGP_TEST_F(NgpHowTo, ReuseNgpField)
{
  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  int numStates = 1;
  stk::mesh::Field<short>                  &shortStkField = create_field_with_num_states_and_init(get_meta(), "field01", numStates, (short)0);
  stk::mesh::Field<unsigned short>        &ushortStkField = create_field_with_num_states_and_init(get_meta(), "field02", numStates, (unsigned short)0);
  stk::mesh::Field<int>                      &intStkField = create_field_with_num_states_and_init(get_meta(), "field03", numStates, (int)0);
  stk::mesh::Field<unsigned int>            &uintStkField = create_field_with_num_states_and_init(get_meta(), "field04", numStates, (unsigned int)0);
  stk::mesh::Field<long>                    &longStkField = create_field_with_num_states_and_init(get_meta(), "field05", numStates, (long)0);
  stk::mesh::Field<unsigned long>          &ulongStkField = create_field_with_num_states_and_init(get_meta(), "field06", numStates, (unsigned long)0);
  stk::mesh::Field<long long>           &longLongStkField = create_field_with_num_states_and_init(get_meta(), "field07", numStates, (long long)0);
  stk::mesh::Field<unsigned long long> &ulongLongStkField = create_field_with_num_states_and_init(get_meta(), "field08", numStates, (unsigned long long)0);
  stk::mesh::Field<float>                  &floatStkField = create_field_with_num_states_and_init(get_meta(), "field09", numStates, (float)0);
  stk::mesh::Field<double>                &doubleStkField = create_field_with_num_states_and_init(get_meta(), "field10", numStates, (double)0);

  stk::io::fill_mesh("generated:1x1x4", get_bulk());
  auto fields = get_meta().get_fields();
  for(auto field : fields) {
    field->sync_to_device();
  }

  {
    fill_field_on_device(get_bulk(), shortStkField, (short)42);
    check_field_on_device(get_bulk(), shortStkField, (short)42);

    fill_field_on_device(get_bulk(), ushortStkField, (unsigned short)43);
    check_field_on_device(get_bulk(), ushortStkField, (unsigned short)43);

    fill_field_on_device(get_bulk(), intStkField, (int)44);
    check_field_on_device(get_bulk(), intStkField, (int)44);

    fill_field_on_device(get_bulk(), uintStkField, (unsigned int)45);
    check_field_on_device(get_bulk(), uintStkField, (unsigned int)45);

    fill_field_on_device(get_bulk(), longStkField, (long)46);
    check_field_on_device(get_bulk(), longStkField, (long)46);

    fill_field_on_device(get_bulk(), ulongStkField, (unsigned long)47);
    check_field_on_device(get_bulk(), ulongStkField, (unsigned long)47);

    fill_field_on_device(get_bulk(), longLongStkField, (long long)48);
    check_field_on_device(get_bulk(), longLongStkField, (long long)48);

    fill_field_on_device(get_bulk(), ulongLongStkField, (unsigned long long)49);
    check_field_on_device(get_bulk(), ulongLongStkField, (unsigned long long)49);

    fill_field_on_device(get_bulk(), floatStkField, (float)3.14);
    check_field_on_device(get_bulk(), floatStkField, (float)3.14);

    fill_field_on_device(get_bulk(), doubleStkField, (double)3.141);
    check_field_on_device(get_bulk(), doubleStkField, (double)3.141);
  }

  //now check field on host to see if FieldManager::clear_fields() (called by FieldManager dtor)
  //sync'd device fields back to host.
  //Only do this check if cuda, because if not cuda then host == device.
#ifdef STK_ENABLE_GPU
  check_field_on_host(get_bulk(), doubleStkField, (double)3.141);
#endif
}


void run_part_membership_test(const stk::mesh::BulkData& bulk, stk::mesh::PartOrdinal partOrdinal)
{
  const stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>::member_type TeamHandleType;
  const auto& teamPolicy = stk::ngp::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
                                                                                 Kokkos::AUTO);

  Kokkos::parallel_for(teamPolicy,
                       KOKKOS_LAMBDA(const TeamHandleType& team)
                       {
                         const stk::mesh::NgpMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK,
                         team.league_rank());
                         unsigned numElems = bucket.size();

                         Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems),
                         [&ngpMesh, &bucket, &partOrdinal](const int& i) {
                           stk::mesh::Entity elem = bucket[i];
                           stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
                           stk::mesh::NgpMesh::ConnectedNodes nodes =
                           ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);

                           for (unsigned nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++) {
                             stk::mesh::Entity node = nodes[nodeIndex];
                             stk::mesh::FastMeshIndex meshIndex = ngpMesh.fast_mesh_index(node);
                             const stk::mesh::NgpMesh::BucketType& nodeBucket =
                             ngpMesh.get_bucket(stk::topology::NODE_RANK, meshIndex.bucket_id);
                             if (ngpMesh.identifier(node) == 1 || ngpMesh.identifier(node) == 2) {
                               STK_NGP_ThrowRequire(nodeBucket.member(partOrdinal) == true);
                             } else {
                               STK_NGP_ThrowRequire(nodeBucket.member(partOrdinal) == false);
                             }
                           }
                         });
                       });
}

TEST_F(NgpHowTo, checkPartMembership)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<double>(stk::topology::NODE_RANK, "myField");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);

  stk::mesh::Part& testPart = get_meta().declare_part("testPart", stk::topology::NODE_RANK);

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1u);
  stk::mesh::Entity node2 = get_bulk().get_entity(stk::topology::NODE_RANK, 2u);

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(node1, stk::mesh::PartVector{&testPart}, {});
  get_bulk().change_entity_parts(node2, stk::mesh::PartVector{&testPart}, {});
  get_bulk().modification_end();

  run_part_membership_test(get_bulk(), testPart.mesh_meta_data_ordinal());
}

void fill_ngp_field(const stk::mesh::NgpMesh& ngpMesh, stk::mesh::EntityRank rank, const stk::mesh::MetaData& meta, stk::mesh::NgpField<int>& ngpField, int fieldVal)
{
  //BEGINNgpMeshIndexUsage
  stk::mesh::for_each_entity_run(ngpMesh, rank, meta.universal_part(), KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   ngpField(entity, 0) = fieldVal;
                                 });
  //ENDNgpMeshIndexUsage
}

TEST(NgpMesh, meshIndices)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;
  unsigned numStates = 1;

  const int init = 1;
  stk::mesh::Field<int> &field = meta.declare_field<int>(rank, "field_1", numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &init);

  stk::io::fill_mesh("generated:1x1x1", *bulk);
  field.sync_to_device();
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(*bulk);
  stk::mesh::NgpField<int> & ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  int fieldVal = 5;

  fill_ngp_field(ngpMesh, rank, meta, ngpField, fieldVal);

  ngpField.modify_on_device();
  ngpField.sync_to_host();

  stk::mesh::EntityId id = 1;
  stk::mesh::Entity entity = bulk->get_entity(rank, id);
  int* data = stk::mesh::field_data(field, entity);

  ASSERT_EQ(fieldVal, data[0]);
}

}

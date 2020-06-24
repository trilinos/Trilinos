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
#include "NgpUnitTestUtils.hpp"

namespace {

using IntDualViewType = Kokkos::DualView<int*, stk::mesh::ExecSpace>;

void set_field_on_device_and_copy_back(stk::mesh::BulkData &bulk,
                                       stk::mesh::EntityRank rank,
                                       stk::mesh::Part &quadPart,
                                       stk::mesh::Field<double> &quadField,
                                       double fieldVal)
{
  stk::mesh::NgpField<double>& ngpQuadField = stk::mesh::get_updated_ngp_field<double>(quadField);
  EXPECT_EQ(quadField.mesh_meta_data_ordinal(), ngpQuadField.get_ordinal());

  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  EXPECT_EQ(bulk.mesh_meta_data().spatial_dimension(), ngpMesh.get_spatial_dimension());

  stk::mesh::for_each_entity_run(ngpMesh, rank, quadPart,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   ngpQuadField(entity, 0) = fieldVal;
                                 });
  ngpQuadField.modify_on_device();
  ngpQuadField.sync_to_host();
}

class NgpHowTo : public stk::unit_test_util::MeshFixture 
{
public:
  void setup_test_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
    auto &shellQuadField = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField");
    double init = 0.0;
    stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, &init);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
        0,2,SHELL_QUAD_4,5,6,7,8";
        stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(NgpHowTo, loopOverSubsetOfMesh)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
  auto &shellQuadField = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, &init);
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
      0,2,SHELL_QUAD_4,5,6,7,8";
      stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  double fieldVal = 13.0;
  set_field_on_device_and_copy_back(get_bulk(), stk::topology::ELEM_RANK, shellQuadPart, shellQuadField, fieldVal);

  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, shellQuadPart))
  {
    for(stk::mesh::Entity elem : *bucket)
    {
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(shellQuadField, elem));
    }
  }
}

template<typename MeshType>
void test_mesh_up_to_date(stk::mesh::BulkData& bulk)
{
  //BEGINNgpMeshUpToDate
  MeshType& ngpMesh = bulk.get_updated_ngp_mesh();
  EXPECT_TRUE(ngpMesh.is_up_to_date());

  bulk.modification_begin();
  bulk.modification_end();

  MeshType& newNgpMesh = bulk.get_updated_ngp_mesh();
  EXPECT_TRUE(newNgpMesh.is_up_to_date());
  //ENDNgpMeshUpToDate
}

TEST_F(NgpHowTo, checkIfUpToDate)
{
  setup_test_mesh();
  test_mesh_up_to_date<stk::mesh::NgpMesh>(get_bulk());
}

template <typename FieldType>
void test_field_on_subset_of_mesh(const stk::mesh::BulkData& bulk, const FieldType& field,
                                  stk::mesh::PartOrdinal partThatHasField,
                                  stk::mesh::PartOrdinal partThatDoesntHaveField)
{
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           if (bucket.member(partThatHasField)) {
                             NGP_ThrowRequire(field.get_num_components_per_entity(elemIndex) > 0);
                           }
                           if (bucket.member(partThatDoesntHaveField)) {
                             NGP_ThrowRequire(field.get_num_components_per_entity(elemIndex) == 0);
                           }
                         });
                       });
}

TEST_F(NgpHowTo, fieldOnSubsetOfMesh)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &shellQuadPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4);
  const stk::mesh::Part &hex8Part = get_meta().get_topology_root_part(stk::topology::HEX_8);

  auto &shellQuadField = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(shellQuadField, shellQuadPart, &init);
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
      0,2,SHELL_QUAD_4,5,6,7,8";
      stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  double fieldVal = 13.0;
  set_field_on_device_and_copy_back(get_bulk(), stk::topology::ELEM_RANK, shellQuadPart, shellQuadField, fieldVal);

  stk::mesh::NgpField<double> & ngpShellFieldAdapter = stk::mesh::get_updated_ngp_field<double>(shellQuadField);

  test_field_on_subset_of_mesh(get_bulk(), ngpShellFieldAdapter,
                               shellQuadPart.mesh_meta_data_ordinal(), hex8Part.mesh_meta_data_ordinal());
}

TEST_F(NgpHowTo, loopOverAllMeshNodes)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  double fieldVal = 13.0;
  set_field_on_device_and_copy_back(get_bulk(), stk::topology::NODE_RANK, get_meta().universal_part(), field, fieldVal);

  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().universal_part()))
    for(stk::mesh::Entity node : *bucket)
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(field, node));
}

TEST_F(NgpHowTo, loopOverMeshFaces)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Part &facePart = get_meta().declare_part("facePart", stk::topology::FACE_RANK);
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::FACE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, facePart, &init);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {&facePart});

  double fieldVal = 13.0;
  set_field_on_device_and_copy_back(get_bulk(), stk::topology::FACE_RANK, facePart, field, fieldVal);

  for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::FACE_RANK, get_meta().universal_part()))
    for(stk::mesh::Entity node : *bucket)
      EXPECT_EQ(fieldVal, *stk::mesh::field_data(field, node));
}

void run_connected_node_test(const stk::mesh::BulkData& bulk)
{
  stk::topology elemTopo = stk::topology::HEX_8;

  stk::mesh::EntityVector elems;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
  EXPECT_EQ(1u, elems.size());
  stk::mesh::Entity node0 = bulk.begin_nodes(elems[0])[0];
  stk::mesh::Entity node7 = bulk.begin_nodes(elems[0])[7];

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);
                           stk::topology bucketTopo = bucket.topology();
                           NGP_ThrowRequire(elemTopo == bucketTopo);
                           NGP_ThrowRequire(nodes.size() == bucketTopo.num_nodes());
                           NGP_ThrowRequire(node0 == nodes[0]);
                           NGP_ThrowRequire(node7 == nodes[7]);

                           stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[0]);
                           stk::mesh::NgpMesh::ConnectedEntities node0_elems = ngpMesh.get_elements(stk::topology::NODE_RANK,
                                                                                                    nodeIndex);
                           NGP_ThrowRequire(1 == node0_elems.size());
                           NGP_ThrowRequire(node0_elems[0] == elem);
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
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
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, bucketCapacity);
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  run_connected_node_test(get_bulk());
}

void run_id_test(const stk::mesh::BulkData& bulk)
{
  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           NGP_ThrowRequire(elemTopo == bucketTopo);
                           NGP_ThrowRequire(faces.size() == bucketTopo.num_faces());
                           NGP_ThrowRequire(face0 == faces[0]);
                           NGP_ThrowRequire(face1 == faces[1]);
                           NGP_ThrowRequire(face2 == faces[2]);
                           NGP_ThrowRequire(face3 == faces[3]);
                           NGP_ThrowRequire(face4 == faces[4]);
                           NGP_ThrowRequire(face5 == faces[5]);

                           stk::mesh::NgpMesh::ConnectedEntities edges = ngpMesh.get_edges(stk::topology::ELEM_RANK, elemIndex);
                           NGP_ThrowRequire(0 == edges.size());

                           stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
                           stk::mesh::NgpMesh::ConnectedEntities face0_elems = ngpMesh.get_elements(stk::topology::FACE_RANK,
                                                                                                    faceIndex);
                           NGP_ThrowRequire(1 == face0_elems.size());
                           NGP_ThrowRequire(face0_elems[0] == elem);
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFaces)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  setup_mesh("generated:1x1x1|sideset:xXyYzZ", stk::mesh::BulkData::NO_AUTO_AURA);

  run_connected_face_test(get_bulk());
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

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           NGP_ThrowRequire(elemTopo == bucketTopo);
                           NGP_ThrowRequire(ordinals.size() == bucketTopo.num_faces());
                           NGP_ThrowRequire(ordinal0 == ordinals[0]);
                           NGP_ThrowRequire(ordinal1 == ordinals[1]);
                           NGP_ThrowRequire(ordinal2 == ordinals[2]);
                           NGP_ThrowRequire(ordinal3 == ordinals[3]);
                           NGP_ThrowRequire(ordinal4 == ordinals[4]);
                           NGP_ThrowRequire(ordinal5 == ordinals[5]);

                           stk::mesh::NgpMesh::ConnectedOrdinals edgeOrdinals = ngpMesh.get_edge_ordinals(stk::topology::ELEM_RANK,
                                                                                                          elemIndex);
                           NGP_ThrowRequire(0 == edgeOrdinals.size());

                           stk::mesh::Entity face = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex)[0];
                           stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(face);
                           stk::mesh::NgpMesh::ConnectedOrdinals face0_elemOrdinals = ngpMesh.get_element_ordinals(stk::topology::FACE_RANK, faceIndex);
                           NGP_ThrowRequire(1 == face0_elemOrdinals.size());
                           NGP_ThrowRequire(face0_elemOrdinals[0] == ordinal0);
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFaceOrdinals)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  setup_mesh("generated:1x1x1|sideset:xXyYzZ", stk::mesh::BulkData::NO_AUTO_AURA);

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

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           NGP_ThrowRequire(elemTopo == bucketTopo);
                           NGP_ThrowRequire(facePermutations.size() == bucketTopo.num_faces());
                           NGP_ThrowRequire(edgePermutations.size() == bucketTopo.num_edges());

                           for(unsigned j=0; j<numFacesPerHex; ++j) {
                             NGP_ThrowRequire(p.elemFacePermutations[i][j] == facePermutations[j]);
                           }
                           for(unsigned j=0; j<numEdgesPerHex; ++j) {
                             NGP_ThrowRequire(p.elemEdgePermutations[i][j] == edgePermutations[j]);
                           }
                           if (i == 0) {
                             stk::mesh::Entity face = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex)[0];
                             stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(face);
                             stk::mesh::NgpMesh::Permutations face0_elemPermutations = ngpMesh.get_element_permutations(stk::topology::FACE_RANK, faceIndex);
                             NGP_ThrowRequire(1 == face0_elemPermutations.size());
                             NGP_ThrowRequire(face0_elemPermutations[0] == p.elemFacePermutations[0][0]);
                           }
                         });
                       });
}

TEST_F(NgpHowTo, loopOverElemFacePermutations)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  stk::mesh::Part& sidePart = get_meta().declare_part("SidePart");
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::NO_AUTO_AURA);
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
  IntDualViewType result = ngp_unit_test_utils::create_dualview<IntDualViewType>("result",numResults);
  enum {ELEM_FACE_CHECK = 0, FACE_NODE_CHECK = 1};

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           NGP_ThrowRequire(elemTopo == bucketTopo);
                           int numFaces = faces.size();
                           stk::mesh::atomic_add(&result.d_view(ELEM_FACE_CHECK), numFaces);

                           if (numFaces == 1)
                           {
                             stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
                             stk::mesh::NgpMesh::ConnectedEntities faceNodes = ngpMesh.get_connected_entities(stk::topology::FACE_RANK, faceIndex, stk::topology::NODE_RANK);

                             unsigned faceOrdinal = 0; //when we add ordinals we can fix this. But for a hex all faces have the same topology anyway...
                             stk::topology faceTopo = elemTopo.face_topology(faceOrdinal);

                             NGP_ThrowRequire(faceNodes.size() == faceTopo.num_nodes());
                             stk::mesh::atomic_add(&result.d_view(FACE_NODE_CHECK), 1);
                           }
                         });
                       });

  result.modify<IntDualViewType::execution_space>();
  result.sync<IntDualViewType::host_mirror_space>();

  EXPECT_EQ(2, result.h_view(ELEM_FACE_CHECK)); //expected 2 elements that had faces
  EXPECT_EQ(2, result.h_view(FACE_NODE_CHECK)); //expected 2 faces that had nodes
}

TEST_F(NgpHowTo, anotherElemFacesTest)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
  setup_mesh("generated:1x1x4|sideset:zZ", stk::mesh::BulkData::NO_AUTO_AURA);

  run_another_connected_face_test(get_bulk());
}

void test_view_of_fields(const stk::mesh::BulkData& bulk,
                         stk::mesh::Field<double>& field1,
                         stk::mesh::Field<double>& field2)
{
  using FieldViewType = Kokkos::View<stk::mesh::NgpField<double>*,stk::mesh::MemSpace>;

  FieldViewType fields(Kokkos::ViewAllocateWithoutInitializing("fields"),2);
  FieldViewType::HostMirror hostFields = Kokkos::create_mirror_view(fields);

  Kokkos::parallel_for(2,
                       KOKKOS_LAMBDA(const unsigned& i)
                       {
                         new (&fields(i)) stk::mesh::NgpField<double>();
                       });

  hostFields(0) = stk::mesh::NgpField<double>(bulk, field1);
  hostFields(1) = stk::mesh::NgpField<double>(bulk, field2);

  Kokkos::deep_copy(fields, hostFields);

  unsigned numResults = 2;
  IntDualViewType result = ngp_unit_test_utils::create_dualview<IntDualViewType>("result",numResults);

  Kokkos::parallel_for(2,
                       KOKKOS_LAMBDA(const unsigned& i)
                       {
                         result.d_view(i) = fields(i).get_ordinal() == i ? 1 : 0;
                       });

  result.modify<IntDualViewType::execution_space>();
  result.sync<IntDualViewType::host_mirror_space>();

  EXPECT_EQ(1, result.h_view(0));
  EXPECT_EQ(1, result.h_view(1));
}

// Disabled because stk::mesh::NgpField now contains a Kokkos::DualView, which is not
// properly constructible on the device in the old version of Kokkos that we
// currently have in Sierra.  Versions after at least 2018-12-10 work fine, so
// this can be re-enabled after our next Trilinos pull.
TEST_F(NgpHowTo, DISABLED_viewOfFields)
{
  auto &field1 = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField1");
  auto &field2 = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField2");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field1, get_meta().universal_part(), &init);
  stk::mesh::put_field_on_mesh(field2, get_meta().universal_part(), &init);
  setup_mesh("generated:1x1x4|sideset:zZ", stk::mesh::BulkData::NO_AUTO_AURA);

  test_view_of_fields(get_bulk(), field1, field2);
}

void test_ngp_mesh_construction(const stk::mesh::BulkData& bulk)
{
  size_t numHostElemBuckets = bulk.buckets(stk::topology::ELEM_RANK).size();
  std::vector<size_t> counts(bulk.mesh_meta_data().entity_rank_count(), 0);
  stk::mesh::count_entities(bulk.mesh_meta_data().locally_owned_part(), bulk, counts);
  size_t numElements = counts[stk::topology::ELEM_RANK];

  double startTime = stk::wall_time();

  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

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

  stk::mesh::create_exposed_block_boundary_sides(get_bulk(), meta.universal_part(), {&boundaryPart});
  stk::mesh::create_interior_block_boundary_sides(get_bulk(), meta.universal_part(), {&boundaryPart});

  test_ngp_mesh_construction(get_bulk());
}

unsigned count_num_elems(stk::mesh::NgpMesh ngpMesh,
                         stk::mesh::NgpField<int> ngpField,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Part &part)
{
  Kokkos::View<unsigned *, stk::mesh::MemSpace> numElems("numElems", 1);
  stk::mesh::for_each_entity_run(ngpMesh, rank, part,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   unsigned fieldValue = static_cast<unsigned>(ngpField(entity, 0));
                                   Kokkos::atomic_add(&numElems(0), fieldValue);
                                 });
  Kokkos::View<unsigned *, stk::mesh::MemSpace>::HostMirror numElemsHost =
      Kokkos::create_mirror_view(numElems);
  Kokkos::deep_copy(numElemsHost, numElems);
  return numElemsHost(0);
}

void set_num_elems_in_field_on_device_and_copy_back(stk::mesh::BulkData &bulk,
                                                    stk::mesh::Part &part,
                                                    stk::mesh::Field<int> &field)
{
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
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
  auto &field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, "myField");
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
  auto &field = meta.declare_field<stk::mesh::Field<DataType>>(stk::topology::ELEM_RANK, fieldName, numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &init);
  return field;
}

template <typename DataType>
stk::mesh::Field<DataType, stk::mesh::Cartesian> &create_vector_field_with_num_states_and_init(stk::mesh::MetaData &meta,
                                                                                               const std::string & fieldName,
                                                                                               int numStates,
                                                                                               int fieldDimension,
                                                                                               DataType* init)
{
  auto &field = meta.declare_field<stk::mesh::Field<DataType, stk::mesh::Cartesian>>(stk::topology::ELEM_RANK, fieldName, numStates);
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
  stk::mesh::Field<double> &stkField = create_field_with_num_states_and_init<double>(get_meta(), "myField", numStates, initialValue);
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();

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
  auto& stkField = create_vector_field_with_num_states_and_init<double>(get_meta(),
                                                                        "myField",
                                                                        numStates,
                                                                        fieldDimension,
                                                                        initialValue);
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  const stk::mesh::NgpMesh& ngpMesh = get_bulk().get_updated_ngp_mesh();

  double fieldVal = 1.0;
  ngpField.set_all(ngpMesh, fieldVal);

  double expectedSum = fieldDimension*stk::mesh::count_selected_entities(get_meta().locally_owned_part(),
                                                                         get_bulk().buckets(stk::topology::ELEM_RANK));
  double sum = stk::mesh::get_field_sum(ngpMesh, ngpField, get_meta().locally_owned_part());
  EXPECT_NEAR(expectedSum, sum, 1e-14);
}

void test_vector_field_size_and_values(stk::mesh::BulkData& bulk,
                              const stk::mesh::FieldBase& stkField)
{
  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(stkField);
  const stk::mesh::NgpMesh& ngpMesh = bulk.get_updated_ngp_mesh();

  const unsigned fieldLength = stkField.max_size(stkField.entity_rank());
  const double tol = 1.e-12;
  const stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part("block_1");
  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), block1,
               KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& meshIdx)
               {
                  NGP_EXPECT_NEAR(1.0, ngpField(meshIdx, 0), tol);
                  NGP_EXPECT_NEAR(2.0, ngpField(meshIdx, 1), tol);
                  NGP_EXPECT_NEAR(3.0, ngpField(meshIdx, 2), tol);

                  stk::mesh::EntityFieldData<double> vals = ngpField(meshIdx);
                  NGP_EXPECT_EQ(vals.size(), ngpField.get_num_components_per_entity(meshIdx));
                  NGP_EXPECT_EQ(fieldLength, vals.size());
                  NGP_EXPECT_NEAR(1.0, vals[0], tol);
                  NGP_EXPECT_NEAR(2.0, vals[1], tol);
                  NGP_EXPECT_NEAR(3.0, vals[2], tol);

                  vals[0] += 2.0;
                  vals[1] += 3.0;
                  vals[2] += 4.0;
               });

  stk::mesh::for_each_entity_run(ngpMesh, ngpField.get_rank(), block1,
               KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& meshIdx)
               {
                  const stk::mesh::EntityFieldData<double> vals = ngpField(meshIdx);

                  const double& vals0 = vals[0];
                  const double& vals1 = vals[1];
                  const double& vals2 = vals[2];

                  NGP_EXPECT_NEAR(3.0, vals0, tol);
                  NGP_EXPECT_NEAR(5.0, vals1, tol);
                  NGP_EXPECT_NEAR(7.0, vals2, tol);
               });
}

NGP_TEST_F(NgpHowTo, accessVectorFieldValues)
{
  int numStates = 1;
  double initialValue[] = {1.0, 2.0, 3.0};
  int fieldDimension = 3;
  auto& stkField = create_vector_field_with_num_states_and_init<double>(get_meta(),
                                                                        "myField",
                                                                        numStates,
                                                                        fieldDimension,
                                                                        initialValue);
  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

  test_vector_field_size_and_values(get_bulk(), stkField);
}


class NgpReduceHowTo : public stk::unit_test_util::MeshFixture
{
protected:
  NgpReduceHowTo()
  {
    elemField = &create_field(get_meta());
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
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

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  const int lowestElemID = 1;
  int expectedMinVal = lowestElemID;
  int actualMinVal = stk::mesh::get_field_min(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedMinVal, actualMinVal);
}

TEST_F(NgpReduceHowTo, getMaxFieldValue)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  const int highestElemID = 4;
  int expectedMaxVal = highestElemID;
  int actualMaxVal = stk::mesh::get_field_max(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedMaxVal, actualMaxVal);
}

TEST_F(NgpReduceHowTo, getSumFieldValue)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int numElems = get_num_elems();
  int expectedSum = numElems*(numElems+1)/2;
  int sum_val = stk::mesh::get_field_sum(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part());
  EXPECT_EQ(expectedSum, sum_val);
}
TEST_F(NgpReduceHowTo, minMaxPairWiseReduction)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
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

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
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

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int expectedMin = 1;
  stk::mesh::EntityId expectedMinLoc = 1;
  int expectedMax = get_num_elems();
  stk::mesh::EntityId expectedMaxLoc = 4;
  Kokkos::MinMaxLocScalar<int,stk::mesh::EntityId> minMaxLocVal;
  Kokkos::MinMaxLoc<int,stk::mesh::EntityId> minMaxLoc (minMaxLocVal);
  stk::mesh::get_field_reduction (ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part(), minMaxLoc);
  EXPECT_EQ(expectedMin, minMaxLocVal.min_val);
  EXPECT_EQ(expectedMinLoc, minMaxLocVal.min_loc);
  EXPECT_EQ(expectedMax, minMaxLocVal.max_val);
  EXPECT_EQ(expectedMaxLoc, minMaxLocVal.max_loc);
}
TEST_F(NgpReduceHowTo, minMaxLocReductionThroughAccessor)
{
  if (get_bulk().parallel_size() > 1) return;

  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpElemField = stk::mesh::get_updated_ngp_field<int>(*elemField);
  int expectedMin = 1;
  stk::mesh::EntityId expectedMinLoc = 1;
  int expectedMax = get_num_elems();
  stk::mesh::EntityId expectedMaxLoc = 4;
  Kokkos::MinMaxLocScalar<int,stk::mesh::EntityId> minMaxLocVal;
  Kokkos::MinMaxLoc<int,stk::mesh::EntityId> minMaxLoc (minMaxLocVal);
  stk::mesh::FieldAccessFunctor<stk::mesh::NgpMesh, stk::mesh::NgpField<int>,
                                decltype(minMaxLoc), stk::mesh::identity<int>>accessor(ngpElemField, minMaxLoc);
  stk::mesh::get_field_reduction(ngpMesh, get_bulk().mesh_meta_data().universal_part(), accessor);
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
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
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
  stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::NgpField<T> & ngpField = stk::mesh::get_updated_ngp_field<T>(stkField);
  ngpField.sync_to_device();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK,
                                 bulk.mesh_meta_data().locally_owned_part(),
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity)
                                 {
                                   NGP_EXPECT_EQ(ngpField(entity, 0), expectedFieldValue);
                                 });
}

template <typename T>
void check_field_on_host(const stk::mesh::BulkData & bulk,
                         stk::mesh::Field<T> & stkField,
                         T expectedFieldValue)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stkField.sync_to_host();

  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  for(const stk::mesh::Bucket* bptr : buckets) {
    for(stk::mesh::Entity elem : *bptr) {
      const double* fieldData = stk::mesh::field_data(stkField, elem);
      EXPECT_EQ(*fieldData, expectedFieldValue);
    }
  }
}

NGP_TEST_F(NgpHowTo, ReuseNgpField)
{
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

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

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
#ifdef KOKKOS_ENABLE_CUDA
  check_field_on_host(get_bulk(), doubleStkField, (double)3.141);
#endif
}


void run_part_membership_test(const stk::mesh::BulkData& bulk, stk::mesh::PartOrdinal partOrdinal)
{
  const stk::mesh::NgpMesh & ngpMesh = bulk.get_updated_ngp_mesh();

  typedef Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace, stk::mesh::ScheduleType>::member_type TeamHandleType;
  const auto& teamPolicy = Kokkos::TeamPolicy<stk::mesh::NgpMesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK),
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
                           stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);

                           for(unsigned nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++) {
                             stk::mesh::Entity node = nodes[nodeIndex];
                             stk::mesh::FastMeshIndex meshIndex = ngpMesh.fast_mesh_index(node);
                             const stk::mesh::NgpMesh::BucketType& nodeBucket = ngpMesh.get_bucket(stk::topology::NODE_RANK,
                                                                                                   meshIndex.bucket_id);
                             if(ngpMesh.identifier(node) == 1 || ngpMesh.identifier(node) == 2) {
                               NGP_ThrowRequire(nodeBucket.member(partOrdinal) == true);
                             } else {
                               NGP_ThrowRequire(nodeBucket.member(partOrdinal) == false);
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
  auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
  double init = 0.0;
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);

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

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;
  unsigned numStates = 1;

  const int init = 1;
  stk::mesh::Field<int> &field = meta.declare_field<stk::mesh::Field<int>>(rank, "field_1", numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), &init);

  stk::io::fill_mesh("generated:1x1x1", bulk);
  stk::mesh::NgpMesh& ngpMesh = bulk.get_updated_ngp_mesh();
  stk::mesh::NgpField<int> & ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  int fieldVal = 5;

  fill_ngp_field(ngpMesh, rank, meta, ngpField, fieldVal);

  ngpField.modify_on_device();
  ngpField.sync_to_host();

  stk::mesh::EntityId id = 1;
  stk::mesh::Entity entity = bulk.get_entity(rank, id);
  int* data = stk::mesh::field_data(field, entity);

  ASSERT_EQ(fieldVal, data[0]);
}

//==============================================================================
class FakeEntity {
public:
  STK_FUNCTION
  FakeEntity()
    : m_value(0)
  {
    printf("  FakeEntity: (%lu) Calling default constructor\n", m_value);
  }

  STK_FUNCTION
  explicit FakeEntity(size_t value)
    : m_value(value)
  {
    printf("  FakeEntity: (%lu) Calling constructor\n", m_value);
  }

  STK_FUNCTION
  ~FakeEntity() {
    printf("  FakeEntity: (%lu) Calling destructor\n", m_value);
  }

  STK_FUNCTION
  FakeEntity(const FakeEntity& rhs) {
    printf("  FakeEntity: (%lu) Calling copy constructor\n", rhs.m_value);
    m_value = rhs.m_value;
  }

  STK_FUNCTION
  size_t value() const { return m_value; }

private:
  size_t m_value;
};

using FakeEntityType = Kokkos::View<FakeEntity*, stk::mesh::MemSpace>;

class FakeBucket {
public:
  STK_FUNCTION
  FakeBucket()
    : m_value(0)
  {
    printf("FakeBucket: (%lu) Calling default constructor\n", m_value);
  }

  STK_FUNCTION
  explicit FakeBucket(size_t value)
    : m_value(value)
  {
    printf("FakeBucket: (%lu) Calling constructor\n", m_value);
  }

  STK_FUNCTION
  ~FakeBucket() {
    printf("FakeBucket: (%lu) Calling destructor\n", m_value);
  }

  STK_FUNCTION
  FakeBucket(const FakeBucket& rhs) {
    printf("FakeBucket: (%lu) Calling copy constructor\n", rhs.m_value);
    m_value = rhs.m_value;
    m_innerView = rhs.m_innerView;
  }

  void initialize(size_t numValues) {
    m_innerView = FakeEntityType("Data", numValues);
  }

  STK_FUNCTION
  size_t value(size_t i) const { return m_innerView[i].value(); }

private:
  size_t m_value;
  FakeEntityType m_innerView;
};

using FakeBuckets = Kokkos::View<FakeBucket*, stk::mesh::UVMMemSpace>;

class FakeMesh
{
public:
  STK_FUNCTION
  FakeMesh()
    : m_isInitialized(false),
      m_numBuckets(1),
      m_numEntities(2)
  {
    printf("FakeMesh: Calling default constructor\n");
    update(0);
  }

  STK_FUNCTION
  ~FakeMesh() {
    printf("FakeMesh: Calling destructor\n");
    if (m_fakeBuckets.use_count() == 1) {
      clear();
    }
  }

  STK_FUNCTION
  FakeMesh(const FakeMesh & rhs) {
    printf("FakeMesh: Calling copy constructor\n");
    m_fakeBuckets = rhs.m_fakeBuckets;
    m_isInitialized = rhs.m_isInitialized;
    m_numBuckets = rhs.m_numBuckets;
    m_numEntities = rhs.m_numEntities;
  }

  FakeMesh & operator=(const FakeMesh & rhs) = delete;

  void clear() {
    for (size_t i = 0; i < m_fakeBuckets.size(); ++i) {
      m_fakeBuckets[i].~FakeBucket();
    }
  }

  void fill(size_t iter) {
    m_fakeBuckets = FakeBuckets(Kokkos::ViewAllocateWithoutInitializing("Outer"), m_numBuckets);
    for (size_t i = 0; i < m_numBuckets; ++i) {
      printf("\nFilling buckets: bucket = %lu (iter = %lu)\n", i+1, iter);
      new (&m_fakeBuckets[i]) FakeBucket(i+1);
      m_fakeBuckets[i].initialize(m_numEntities);
    }

  }

  void update(size_t iter) {
    if (m_isInitialized) {
      clear();
    }
    fill(iter);
    m_isInitialized = true;
  }

  STK_FUNCTION
  void do_stuff() const {
    for (size_t i = 0; i < m_numBuckets; ++i) {
      for (size_t j = 0; j < m_numEntities; ++j) {
        printf("Doing stuff: FakeBucket value = %lu\n", m_fakeBuckets[i].value(j));
      }
    }
  }

private:
  FakeBuckets m_fakeBuckets;
  bool m_isInitialized;
  size_t m_numBuckets;
  size_t m_numEntities;
};

void check_fake_mesh_on_device()
{
  FakeMesh fakeMesh;

  const int numMeshMods = 3;
  for (int i = 0; i < numMeshMods; ++i) {
    fakeMesh.update(i+1);
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int idx) {
                           fakeMesh.do_stuff();
                         });
  }

}

TEST(NgpExperiment, DISABLED_FakeMesh)
{
  check_fake_mesh_on_device();
}

//==============================================================================

}

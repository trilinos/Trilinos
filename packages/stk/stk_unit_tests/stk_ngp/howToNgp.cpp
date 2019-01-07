#include <gtest/gtest.h>
#include <stk_ngp/Ngp.hpp>
#include <stk_ngp/NgpAtomics.hpp>
#include <stk_ngp/NgpMultistateField.hpp>
#include <stk_ngp/NgpFieldManager.hpp>
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
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/stk_config.h>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "NgpUnitTestUtils.hpp"

typedef Kokkos::DualView<int*, Kokkos::LayoutRight, ngp::ExecSpace> IntViewType;

extern int gl_argc;
extern char** gl_argv;

void set_field_on_device_and_copy_back(stk::mesh::BulkData &bulk,
                                       stk::mesh::EntityRank rank,
                                       stk::mesh::Part &quadPart,
                                       stk::mesh::Field<double> &quadField,
                                       double fieldVal)
{
    ngp::Field<double> ngpQuadField(bulk, quadField);
    EXPECT_EQ(quadField.mesh_meta_data_ordinal(), ngpQuadField.get_ordinal());

    ngp::Mesh ngpMesh(bulk);
    EXPECT_EQ(bulk.mesh_meta_data().spatial_dimension(), ngpMesh.get_spatial_dimension());

    ngp::for_each_entity_run(ngpMesh, rank, quadPart, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpQuadField.get(entity, 0) = fieldVal;
    });
    ngpQuadField.copy_device_to_host(bulk, quadField);
}

class NgpHowTo : public stk::unit_test_util::MeshFixture {};

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
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    double fieldVal = 13.0;
    set_field_on_device_and_copy_back(get_bulk(), stk::topology::ELEM_RANK, shellQuadPart, shellQuadField, fieldVal);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::ELEM_RANK, shellQuadPart))
        for(stk::mesh::Entity elem : *bucket)
            EXPECT_EQ(fieldVal, *stk::mesh::field_data(shellQuadField, elem));
}

TEST_F(NgpHowTo, loopOverAllMeshNodes)
{
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    auto &field = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "myField");
    double init = 0.0;
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

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
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {&facePart});

    double fieldVal = 13.0;
    set_field_on_device_and_copy_back(get_bulk(), stk::topology::FACE_RANK, facePart, field, fieldVal);

    for(const stk::mesh::Bucket *bucket : get_bulk().get_buckets(stk::topology::FACE_RANK, get_meta().universal_part()))
        for(stk::mesh::Entity node : *bucket)
            EXPECT_EQ(fieldVal, *stk::mesh::field_data(field, node));
}

template<typename MeshType>
void run_connected_node_test(const stk::mesh::BulkData& bulk)
{
    stk::topology elemTopo = stk::topology::HEX_8;

    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
    EXPECT_EQ(1u, elems.size());
    stk::mesh::Entity node0 = bulk.begin_nodes(elems[0])[0];
    stk::mesh::Entity node7 = bulk.begin_nodes(elems[0])[7];

    MeshType ngpMesh(bulk);

    typedef Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK), Kokkos::AUTO);

    Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const typename MeshType::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, team.league_rank());
        unsigned numElems = bucket.size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
        {
            stk::mesh::Entity elem = bucket[i];
            stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
            typename MeshType::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);
            stk::topology bucketTopo = bucket.topology();
            NGP_ThrowRequire(elemTopo == bucketTopo);
            NGP_ThrowRequire(nodes.size() == bucketTopo.num_nodes());
            NGP_ThrowRequire(node0 == nodes[0]);
            NGP_ThrowRequire(node7 == nodes[7]);

            stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[0]);
            typename MeshType::ConnectedEntities node0_elems = ngpMesh.get_elements(stk::topology::NODE_RANK, nodeIndex);
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
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

#ifndef KOKKOS_HAVE_CUDA
    run_connected_node_test<ngp::StkMeshAdapter>(get_bulk());
#endif
    run_connected_node_test<ngp::StaticMesh>(get_bulk());
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

    ngp::StaticMesh ngpMesh(bulk);

    typedef Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK), Kokkos::AUTO);

    Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const ngp::StaticMesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, team.league_rank());
        unsigned numElems = bucket.size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
        {
            stk::mesh::Entity elem = bucket[i];
            stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
            ngp::StaticMesh::ConnectedEntities faces = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex);
            stk::topology bucketTopo = bucket.topology();
            NGP_ThrowRequire(elemTopo == bucketTopo);
            NGP_ThrowRequire(faces.size() == bucketTopo.num_faces());
            NGP_ThrowRequire(face0 == faces[0]);
            NGP_ThrowRequire(face1 == faces[1]);
            NGP_ThrowRequire(face2 == faces[2]);
            NGP_ThrowRequire(face3 == faces[3]);
            NGP_ThrowRequire(face4 == faces[4]);
            NGP_ThrowRequire(face5 == faces[5]);

            ngp::StaticMesh::ConnectedEntities edges = ngpMesh.get_edges(stk::topology::ELEM_RANK, elemIndex);
            NGP_ThrowRequire(0 == edges.size());

            stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
            ngp::StaticMesh::ConnectedEntities face0_elems = ngpMesh.get_elements(stk::topology::FACE_RANK, faceIndex);
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

template<typename MeshType>
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

    MeshType ngpMesh(bulk);

    typedef Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK), Kokkos::AUTO);

    Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const typename MeshType::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, team.league_rank());
        unsigned numElems = bucket.size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
        {
            stk::mesh::Entity elem = bucket[i];
            stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
            typename MeshType::ConnectedOrdinals ordinals = ngpMesh.get_face_ordinals(stk::topology::ELEM_RANK, elemIndex);
            stk::topology bucketTopo = bucket.topology();
            NGP_ThrowRequire(elemTopo == bucketTopo);
            NGP_ThrowRequire(ordinals.size() == bucketTopo.num_faces());
            NGP_ThrowRequire(ordinal0 == ordinals[0]);
            NGP_ThrowRequire(ordinal1 == ordinals[1]);
            NGP_ThrowRequire(ordinal2 == ordinals[2]);
            NGP_ThrowRequire(ordinal3 == ordinals[3]);
            NGP_ThrowRequire(ordinal4 == ordinals[4]);
            NGP_ThrowRequire(ordinal5 == ordinals[5]);

            typename MeshType::ConnectedOrdinals edgeOrdinals = ngpMesh.get_edge_ordinals(stk::topology::ELEM_RANK, elemIndex);
            NGP_ThrowRequire(0 == edgeOrdinals.size());

            stk::mesh::Entity face = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex)[0];
            stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(face);
            typename MeshType::ConnectedOrdinals face0_elemOrdinals = ngpMesh.get_element_ordinals(stk::topology::FACE_RANK, faceIndex);
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

#ifndef KOKKOS_HAVE_CUDA
    run_connected_face_ordinal_test<ngp::StkMeshAdapter>(get_bulk());
#endif
    run_connected_face_ordinal_test<ngp::StaticMesh>(get_bulk());
}

constexpr unsigned numFacesPerHex = 6;
constexpr unsigned numEdgesPerHex = 12;
constexpr unsigned numHexes = 2;

struct PermutationData {
  stk::mesh::Permutation elemFacePermutations[numHexes][numFacesPerHex];
  stk::mesh::Permutation elemEdgePermutations[numHexes][numEdgesPerHex];
};

template<typename MeshType>
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

    MeshType ngpMesh(bulk);

    typedef Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK), Kokkos::AUTO);

    Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const typename MeshType::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, team.league_rank());
        unsigned numElems = bucket.size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
        {
            stk::mesh::Entity elem = bucket[i];
            stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
            typename MeshType::Permutations facePermutations = ngpMesh.get_face_permutations(stk::topology::ELEM_RANK, elemIndex);
            typename MeshType::Permutations edgePermutations = ngpMesh.get_edge_permutations(stk::topology::ELEM_RANK, elemIndex);
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
                typename MeshType::Permutations face0_elemPermutations = ngpMesh.get_element_permutations(stk::topology::FACE_RANK, faceIndex);
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

#ifndef KOKKOS_HAVE_CUDA
    run_connected_face_permutation_test<ngp::StkMeshAdapter>(get_bulk());
#endif
    run_connected_face_permutation_test<ngp::StaticMesh>(get_bulk());
}

void run_another_connected_face_test(const stk::mesh::BulkData& bulk)
{
    stk::topology elemTopo = stk::topology::HEX_8;

    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);
    EXPECT_EQ(4u, elems.size());

    unsigned numResults = 2;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);
    enum {ELEM_FACE_CHECK = 0, FACE_NODE_CHECK = 1};

    ngp::Mesh ngpMesh(bulk);

    typedef Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    const auto& teamPolicy = Kokkos::TeamPolicy<ngp::Mesh::MeshExecSpace>(ngpMesh.num_buckets(stk::topology::ELEM_RANK), Kokkos::AUTO);
    Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const ngp::Mesh::BucketType& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, team.league_rank());
        unsigned numElems = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElems), [&] (const int& i)
        {
            stk::mesh::Entity elem = bucket[i];
            stk::mesh::FastMeshIndex elemIndex = ngpMesh.fast_mesh_index(elem);
            ngp::Mesh::ConnectedEntities faces = ngpMesh.get_faces(stk::topology::ELEM_RANK, elemIndex);
            stk::topology bucketTopo = bucket.topology();
            NGP_ThrowRequire(elemTopo == bucketTopo);
            int numFaces = faces.size();
            ngp::atomic_add(&result.d_view(ELEM_FACE_CHECK), numFaces);

            if (numFaces == 1)
            {
                stk::mesh::FastMeshIndex faceIndex = ngpMesh.fast_mesh_index(faces[0]);
                ngp::Mesh::ConnectedEntities faceNodes = ngpMesh.get_connected_entities(stk::topology::FACE_RANK, faceIndex, stk::topology::NODE_RANK);

                unsigned faceOrdinal = 0; //when we add ordinals we can fix this. But for a hex all faces have the same topology anyway...
                stk::topology faceTopo = elemTopo.face_topology(faceOrdinal);

                NGP_ThrowRequire(faceNodes.size() == faceTopo.num_nodes());
                ngp::atomic_add(&result.d_view(FACE_NODE_CHECK), 1);
            }
        });
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();

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
#ifdef KOKKOS_HAVE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif
  using FieldViewType = Kokkos::View<ngp::Field<double>*,DeviceSpace>;

  FieldViewType fields("fields",2);
  FieldViewType::HostMirror hostFields = Kokkos::create_mirror_view(fields);
  hostFields(0) = ngp::Field<double>(bulk, field1);
  hostFields(1) = ngp::Field<double>(bulk, field2);
  Kokkos::deep_copy(fields, hostFields);

  unsigned numResults = 2;
  IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);
  
  Kokkos::parallel_for(2, KOKKOS_LAMBDA(const unsigned& i)
  {
    result.d_view(i) = fields(i).get_ordinal() == i ? 1 : 0;
  });

  result.modify<IntViewType::execution_space>();
  result.sync<IntViewType::host_mirror_space>();

  EXPECT_EQ(1, result.h_view(0));
  EXPECT_EQ(1, result.h_view(1));
}

TEST_F(NgpHowTo, viewOfFields)
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

  ngp::StaticMesh ngpMesh(bulk);

  double elapsedTime = stk::wall_time() - startTime;
  std::cout << "Time to construct ngp::StaticMesh with "<<numElements<<" elements: "<<elapsedTime << std::endl;

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

unsigned count_num_elems(ngp::Mesh ngpMesh,
                         ngp::Field<int> ngpField,
                         stk::mesh::EntityRank rank,
                         stk::mesh::Part &part)
{
#ifdef KOKKOS_HAVE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif
  Kokkos::View<unsigned *, DeviceSpace> numElems("numElems", 1);
    ngp::for_each_entity_run(ngpMesh, rank, part, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        unsigned fieldValue = static_cast<unsigned>(ngpField.get(entity, 0));
        Kokkos::atomic_add(&numElems(0), fieldValue);
    });
    Kokkos::View<unsigned *, DeviceSpace>::HostMirror numElemsHost =
        Kokkos::create_mirror_view(numElems);
    Kokkos::deep_copy(numElemsHost, numElems);
    return numElemsHost(0);
}

void set_num_elems_in_field_on_device_and_copy_back(stk::mesh::BulkData &bulk,
                                                    stk::mesh::Part &part,
                                                    stk::mesh::Field<int> &field)
{
    ngp::Field<int> ngpField(bulk, field);
    ngp::Mesh ngpMesh(bulk);
    unsigned numElems = count_num_elems(ngpMesh, ngpField, field.entity_rank(), part);
    ngp::for_each_entity_run(ngpMesh, field.entity_rank(), part, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpField.get(entity, 0) = numElems;
    });
    ngpField.copy_device_to_host(bulk, field);
}

TEST_F(NgpHowTo, exerciseAura)
{
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    auto &field = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK, "myField");
    int init = 1;
    stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), &init);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         1,2,HEX_8,5,6,7,8,9,10,11,12";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

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

stk::mesh::Field<int> &create_field(stk::mesh::MetaData &meta)
{
    int numStates = 1;
    int initialValue = -1;
    return create_field_with_num_states_and_init(meta, "myField", numStates, initialValue);
}

void verify_state_new_has_value(stk::mesh::BulkData &bulk,
                                stk::mesh::Field<int>& field,
                                ngp::MultistateField<int> &ngpMultistateField,
                                int np1Value)
{
    ngpMultistateField.copy_device_to_host(bulk, field);
    for(const stk::mesh::Bucket* bucket : bulk.buckets(stk::topology::ELEM_RANK))
        for(stk::mesh::Entity elem : *bucket)
            EXPECT_EQ(np1Value, *static_cast<int*>(stk::mesh::field_data(*field.field_state(stk::mesh::StateNP1), elem)));
}

void set_new_as_old_plus_one_on_device(ngp::Mesh &ngpMesh,
                                       stk::mesh::EntityRank rank,
                                       stk::mesh::Selector sel,
                                       ngp::MultistateField<int> &ngpMultistateField)
{
    int component = 0;
    ngp::for_each_entity_run(ngpMesh, rank, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngp::ConstField<int> stateNField = ngpMultistateField.get_field_old_state(stk::mesh::StateN);
        ngp::Field<int> stateNp1Field = ngpMultistateField.get_field_new_state();
        stateNp1Field.get(entity, component) = stateNField.get(entity, component) + 1;
    });
}

TEST_F(NgpHowTo, useMultistateFields)
{
    int numStates = 2;
    int initialValue = 0;
    stk::mesh::Field<int> &stkField = create_field_with_num_states_and_init(get_meta(), "myField", numStates, initialValue);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

    ngp::MultistateField<int> ngpMultistateField(get_bulk(), stkField);
    ngp::Mesh ngpMesh(get_bulk());

    set_new_as_old_plus_one_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);
    int newStateValue = 1;
    verify_state_new_has_value(get_bulk(), stkField, ngpMultistateField, newStateValue);

    get_bulk().update_field_data_states();
    ngpMultistateField.increment_state();

    set_new_as_old_plus_one_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);
    newStateValue = 2;
    verify_state_new_has_value(get_bulk(), stkField, ngpMultistateField, newStateValue);
}

void set_new_as_old_plus_one_in_convenient_field_on_device(ngp::Mesh &ngpMesh,
                                                           stk::mesh::EntityRank rank,
                                                           stk::mesh::Selector sel,
                                                           ngp::ConvenientMultistateField<int> &ngpMultistateField)
{
    int component = 0;
    ngp::for_each_entity_run(ngpMesh, rank, sel, KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpMultistateField.get_new(entity, component) = ngpMultistateField.get_old(stk::mesh::StateOld, entity, component) + 1;
    });
}

TEST_F(NgpHowTo, useConvenientMultistateFields)
{
    int numStates = 2;
    int initialValue = 0;
    stk::mesh::Field<int> &stkField = create_field_with_num_states_and_init(get_meta(), "myField", numStates, initialValue);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

    ngp::ConvenientMultistateField<int> ngpMultistateField(get_bulk(), stkField);
    ngp::Mesh ngpMesh(get_bulk());

    set_new_as_old_plus_one_in_convenient_field_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);
    int newStateValue = 1;
    verify_state_new_has_value(get_bulk(), stkField, ngpMultistateField, newStateValue);

    get_bulk().update_field_data_states();
    ngpMultistateField.increment_state();

    set_new_as_old_plus_one_in_convenient_field_on_device(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), ngpMultistateField);
    newStateValue = 2;
    verify_state_new_has_value(get_bulk(), stkField, ngpMultistateField, newStateValue);
}

TEST_F(NgpHowTo, setAllFieldValues)
{
    int numStates = 2;
    double initialValue = 0.0;
    stk::mesh::Field<double> &stkField = create_field_with_num_states_and_init<double>(get_meta(), "myField", numStates, initialValue);
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);

    ngp::Field<double> ngpField(get_bulk(), stkField);
    ngp::Mesh ngpMesh(get_bulk());

    double fieldVal = 1.0;
    ngpField.set_all(ngpMesh, fieldVal);

    double expectedSum = 4.0;
    double sum = ngp::get_field_sum(ngpMesh, ngpField, get_meta().universal_part());
    EXPECT_NEAR(expectedSum, sum, 1e-14);
}


class NgpReduceHowTo : public stk::unit_test_util::MeshFixture
{
protected:
    NgpReduceHowTo()
    {
        stk::mesh::Field<int> &elemField = create_field(get_meta());
        setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);
        for(stk::mesh::Entity elem : elems)
        {
            int *fieldData = stk::mesh::field_data(elemField, elem);
            fieldData[0] = get_bulk().identifier(elem);
        }
        ngpMesh = ngp::Mesh(get_bulk());
        ngpElemField = ngp::Field<int>(get_bulk(), elemField);
    }
    int get_num_elems()
    {
        return stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
    }
    ngp::Mesh ngpMesh;
    ngp::Field<int> ngpElemField;
};

int get_min_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_min(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getMinFieldValue)
{
    int expectedMinVal = 1;
    EXPECT_EQ(expectedMinVal, get_min_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}

int get_max_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_max(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getMaxFieldValue)
{
    EXPECT_EQ(get_num_elems(), get_max_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}

int get_sum_field_value(ngp::Mesh &ngpMesh, ngp::Field<int> &ngpElemField, stk::mesh::Selector sel)
{
    return ngp::get_field_sum(ngpMesh, ngpElemField, sel);
}
TEST_F(NgpReduceHowTo, getSumFieldValue)
{
    int numElems = get_num_elems();
    int expectedSum = numElems*(numElems+1)/2;
    EXPECT_EQ(expectedSum, get_sum_field_value(ngpMesh, ngpElemField, get_bulk().mesh_meta_data().universal_part()));
}




template <typename T>
void fill_field_on_device(stk::mesh::BulkData & bulk,
                          ngp::Mesh &mesh,
                          ngp::FieldManager &fieldManager,
                          unsigned fieldOrdinal,
                          T fieldValue)
{
    ngp::Field<T> & field = fieldManager.get_field<T>(fieldOrdinal);

    ngp::for_each_entity_run(mesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        field.get(entity, 0) = fieldValue;
    });
}

template <typename T>
void check_field_on_device(stk::mesh::BulkData & bulk,
                           ngp::Mesh &mesh,
                           ngp::FieldManager &fieldManager,
                           unsigned fieldOrdinal,
                           T expectedFieldValue)
{
    ngp::Field<T> & field = fieldManager.get_field<T>(fieldOrdinal);

    ngp::for_each_entity_run(mesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part(), KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        NGP_ThrowRequire(field.get(entity, 0) == expectedFieldValue);
    });
}

TEST_F(NgpHowTo, ReuseNgpField)
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

    ngp::Mesh ngpMesh(get_bulk());
    ngp::FieldManager fieldManager(get_bulk());
    fieldManager.add_fields( {&shortStkField,
                              &ushortStkField,
                              &intStkField,
                              &uintStkField,
                              &longStkField,
                              &ulongStkField,
                              &longLongStkField,
                              &ulongLongStkField,
                              &floatStkField,
                              &doubleStkField} );

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  shortStkField.mesh_meta_data_ordinal(), (short)42);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, shortStkField.mesh_meta_data_ordinal(), (short)42);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  ushortStkField.mesh_meta_data_ordinal(), (unsigned short)43);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, ushortStkField.mesh_meta_data_ordinal(), (unsigned short)43);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  intStkField.mesh_meta_data_ordinal(), (int)44);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, intStkField.mesh_meta_data_ordinal(), (int)44);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  uintStkField.mesh_meta_data_ordinal(), (unsigned int)45);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, uintStkField.mesh_meta_data_ordinal(), (unsigned int)45);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  longStkField.mesh_meta_data_ordinal(), (long)46);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, longStkField.mesh_meta_data_ordinal(), (long)46);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  ulongStkField.mesh_meta_data_ordinal(), (unsigned long)47);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, ulongStkField.mesh_meta_data_ordinal(), (unsigned long)47);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  longLongStkField.mesh_meta_data_ordinal(), (long long)48);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, longLongStkField.mesh_meta_data_ordinal(), (long long)48);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  ulongLongStkField.mesh_meta_data_ordinal(), (unsigned long long)49);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, ulongLongStkField.mesh_meta_data_ordinal(), (unsigned long long)49);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  floatStkField.mesh_meta_data_ordinal(), (float)3.14);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, floatStkField.mesh_meta_data_ordinal(), (float)3.14);

    fill_field_on_device(get_bulk(), ngpMesh, fieldManager,  doubleStkField.mesh_meta_data_ordinal(), (double)3.141);
    check_field_on_device(get_bulk(), ngpMesh, fieldManager, doubleStkField.mesh_meta_data_ordinal(), (double)3.141);
}

#include <gtest/gtest.h>
#include <stk_ngp/NgpSpaces.hpp>
#include <stk_ngp/Ngp.hpp>
#include <stk_ngp/NgpDynamicMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include "NgpUnitTestUtils.hpp"

typedef Kokkos::DualView<int*, Kokkos::LayoutRight, ngp::ExecSpace> IntViewType;

class NgpMeshModImpl : public stk::unit_test_util::MeshFixture
{
protected:
  void initialize_mesh(const std::string& meshDesc)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
  }
};

void test_bucket_topology(stk::mesh::BulkData& bulk)
{
    ngp::DynamicMesh ngpMesh(bulk);

    unsigned numResults = 6;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);

    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    stk::topology hostTopo = bulk.bucket(elem1).topology();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      stk::mesh::FastMeshIndex fmi = ngpMesh.device_mesh_index(elem1);
      const ngp::DynamicBucket& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, fmi.bucket_id);
      result.d_view(0) = bucket.topology() == hostTopo ? 1 : 0;
      result.d_view(1) = bucket.topology() == stk::topology::HEX_8 ? 1 : 0;
      result.d_view(2) = bucket.topology().num_nodes() == 8 ? 1 : 0;
      result.d_view(3) = bucket.topology().num_sides() == 6 ? 1 : 0;
      result.d_view(4) = bucket.topology().side_rank() == stk::topology::FACE_RANK ? 1 : 0;
      result.d_view(5) = bucket.topology().side_topology() == stk::topology::QUAD_4 ? 1 : 0;
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();

    EXPECT_EQ(numResults, result.h_view.size());

    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, result.h_view(i)) << "failed on result "<<i;
    }
}

void test_bucket_part_info(stk::mesh::BulkData& bulk)
{
    ngp::DynamicMesh ngpMesh(bulk);

    unsigned numResults = 6;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);

    stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, 8);
    unsigned hostBucketId = bulk.bucket(node8).bucket_id();
    size_t hostNumPartsForBucket = bulk.bucket(node8).supersets().size();
    int firstPartOrd = bulk.bucket(node8).supersets()[0]->mesh_meta_data_ordinal();
    int lastPartOrd = bulk.bucket(node8).supersets().back()->mesh_meta_data_ordinal();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      stk::mesh::FastMeshIndex fmi = ngpMesh.device_mesh_index(node8);
      const ngp::DynamicBucket& bucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, fmi.bucket_id);
      result.d_view(0) = bucket.bucket_id()==hostBucketId ? 1 : 0;
      result.d_view(1) = bucket.get_num_parts()==hostNumPartsForBucket ? 1 : 0;
      stk::mesh::Entity deviceNode8 = bucket[fmi.bucket_ord];
      result.d_view(2) = deviceNode8==node8 ? 1 : 0;
      result.d_view(3) = bucket.is_member(firstPartOrd) ? 1 : 0;
      result.d_view(4) = bucket.is_member(lastPartOrd) ? 1 : 0;
      unsigned numElemBuckets = ngpMesh.num_buckets(stk::topology::ELEM_RANK);
      NGP_ThrowRequire(2u == numElemBuckets);
      const ngp::DynamicBucket& bucket0 = ngpMesh.get_bucket(stk::topology::ELEM_RANK, 0);
      const ngp::UnmanagedPartOrdViewType& parts0 = bucket0.get_parts();
      const ngp::DynamicBucket& bucket1 = ngpMesh.get_bucket(stk::topology::ELEM_RANK, 1);
      const ngp::UnmanagedPartOrdViewType& parts1 = bucket1.get_parts();
      unsigned numParts0 = parts0(0);
      unsigned numParts1 = parts1(0);
      NGP_ThrowRequire(numParts0 == numParts1);
      result.d_view(5) = ngp::all_parts_match(parts0, parts1) ? 0 : 1; //expecting false
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();

    EXPECT_EQ(numResults, result.h_view.size());

    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, result.h_view(i)) << "failed on result "<<i;
    }
}

void test_find_bucket_with_parts(stk::mesh::BulkData& bulk)
{
    ngp::DynamicMesh ngpMesh(bulk);

    unsigned numResults = 3;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);

    stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, 8);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      stk::mesh::FastMeshIndex fmi = ngpMesh.device_mesh_index(node8);
      const ngp::DynamicBucket& bucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, fmi.bucket_id);
      const ngp::UnmanagedPartOrdViewType& parts = bucket.get_parts();
      int bktIndex = ngpMesh.find_bucket_with_parts(stk::topology::NODE_RANK, parts);
      result.d_view(0) = bktIndex == (int)fmi.bucket_id ? 1 : 0;
      ngp::UnmanagedPartOrdViewType emptyParts;
      result.d_view(1) = ngp::all_parts_match(emptyParts, parts) ? 0 : 1; //expecting false
      result.d_view(2) = ngp::all_parts_match(emptyParts, emptyParts) ? 1 : 0; //expecting true
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();

    EXPECT_EQ(numResults, result.h_view.size());

    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, result.h_view(i)) << "failed on result "<<i;
    }
}

void test_create_new_part_ord_view(stk::mesh::BulkData& bulk)
{
    ngp::DynamicMesh ngpMesh(bulk);

    unsigned numResults = 1;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);

    unsigned biggestPartOrd = bulk.mesh_meta_data().get_parts().back()->mesh_meta_data_ordinal();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      const ngp::DynamicBucket& bucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, 0);
      unsigned newPartOrd = biggestPartOrd + 1;
      const ngp::UnmanagedPartOrdViewType& oldParts = bucket.get_parts();
      ngp::UnmanagedPartOrdViewType newParts = ngpMesh.get_new_part_ords(newPartOrd, oldParts);
      unsigned newNumParts = newParts(0);
      unsigned oldNumParts = bucket.get_num_parts();
      result.d_view(0) = newNumParts == oldNumParts+1 ? 1 : 0;
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();
    EXPECT_EQ(numResults, result.h_view.size());

    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, result.h_view(i)) << "failed on result "<<i;
    }
}

void test_add_bucket_on_device(stk::mesh::BulkData& bulk)
{
    ngp::DynamicMesh ngpMesh(bulk);

    unsigned numResults = 2;
    IntViewType result = ngp_unit_test_utils::create_dualview<IntViewType>("result",numResults);

    unsigned biggestPartOrd = bulk.mesh_meta_data().get_parts().back()->mesh_meta_data_ordinal();

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      unsigned numNodeBuckets = ngpMesh.num_buckets(stk::topology::NODE_RANK);
      const ngp::DynamicBucket& bucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, 0);
      unsigned newPartOrd = biggestPartOrd + 1;
      const ngp::UnmanagedPartOrdViewType& oldParts = bucket.get_parts();
      ngp::UnmanagedPartOrdViewType newParts = ngpMesh.get_new_part_ords(newPartOrd, oldParts);
      int newBucketIndex = ngpMesh.add_bucket(stk::topology::NODE_RANK, bucket.topology(), newParts);
      int bktIndex = ngpMesh.find_bucket_with_parts(stk::topology::NODE_RANK, newParts);
      result.d_view(0) = ((bktIndex == newBucketIndex) && (bktIndex != -1)) ? 1 : 0;
      unsigned newNumNodeBuckets = ngpMesh.num_buckets(stk::topology::NODE_RANK);
      result.d_view(1) = newNumNodeBuckets == numNodeBuckets+1 ? 1 : 0;
    });

    result.modify<IntViewType::execution_space>();
    result.sync<IntViewType::host_mirror_space>();
    EXPECT_EQ(numResults, result.h_view.size());

    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, result.h_view(i)) << "failed on result "<<i;
    }
}

TEST_F(NgpMeshModImpl, bucketTopology)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    initialize_mesh("0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                    "0,2,SHELL_QUAD_4,5,6,7,8");

    test_bucket_topology(get_bulk());
}

TEST_F(NgpMeshModImpl, bucketPartInfo)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    initialize_mesh("0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                    "0,2,SHELL_QUAD_4,5,6,7,8");

    test_bucket_part_info(get_bulk());
}

TEST_F(NgpMeshModImpl, findBucketWithParts)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    initialize_mesh("0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                    "0,2,SHELL_QUAD_4,5,6,7,8");

    test_find_bucket_with_parts(get_bulk());
}

TEST_F(NgpMeshModImpl, createNewPartOrdView)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    initialize_mesh("0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                    "0,2,SHELL_QUAD_4,5,6,7,8");

    test_create_new_part_ord_view(get_bulk());
}

TEST_F(NgpMeshModImpl, addBucketOnDevice)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    initialize_mesh("0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                    "0,2,SHELL_QUAD_4,5,6,7,8");

    test_add_bucket_on_device(get_bulk());
}


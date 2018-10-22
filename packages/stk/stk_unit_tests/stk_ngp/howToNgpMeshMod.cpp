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

typedef Kokkos::View<int*, Kokkos::LayoutRight, ngp::MemSpace> IntViewType;

class NgpMeshModHowTo : public stk::unit_test_util::MeshFixture {};

void test_change_entity_parts(stk::mesh::BulkData& bulk, unsigned active_ordinal)
{

    ngp::DynamicMesh ngpMesh(bulk);

    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);

    unsigned numResults = 6;
    IntViewType result("result",numResults);
    IntViewType::HostMirror hostResult = Kokkos::create_mirror_view(result);
    for(unsigned i=0; i<numResults; ++i) {
      hostResult(i) = 0;
    }
    Kokkos::deep_copy(result, hostResult);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) {
      stk::mesh::FastMeshIndex oldMeshIndex = ngpMesh.device_mesh_index(elem1);
      const ngp::DynamicBucket& oldBucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, oldMeshIndex.bucket_id);
      unsigned oldBucketSize = oldBucket.size();

      ngpMesh.change_entity_parts(elem1, active_ordinal);

      stk::mesh::FastMeshIndex newMeshIndex = ngpMesh.device_mesh_index(elem1);
      const ngp::DynamicBucket& newBucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, newMeshIndex.bucket_id);
      result(0) = newBucket.is_member(active_ordinal) ? 1 : 0;
      result(1) = newBucket.size()==1 ? 1 : 0;
      result(2) = newBucket[0] == elem1 ? 1 : 0;
      result(3) = oldBucket.size() == oldBucketSize-1 ? 1 : 0;
      result(4) = newBucket.topology() == stk::topology::HEX_8 ? 1 : 0;
      result(5) = oldBucket.topology() == stk::topology::HEX_8 ? 1 : 0;
    });

    Kokkos::deep_copy(hostResult, result);
    EXPECT_EQ(numResults, hostResult.size());
    for(unsigned i=0; i<numResults; ++i) {
      EXPECT_EQ(1, hostResult(i)) << "failed for result "<<i<<std::endl;
    }
}

TEST_F(NgpMeshModHowTo, changeEntityParts)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Part &active = get_meta().declare_part("active", stk::topology::ELEM_RANK);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,SHELL_QUAD_4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());

    test_change_entity_parts(get_bulk(), active.mesh_meta_data_ordinal());
}


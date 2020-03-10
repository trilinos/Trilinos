#include <gtest/gtest.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpAtomics.hpp>
#include <stk_mesh/base/NgpMultistateField.hpp>
#include <stk_mesh/base/NgpFieldManager.hpp>
#include <stk_mesh/base/NgpReductions.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace {

class UpdateNgpMesh : public stk::unit_test_util::MeshFixture
{
public:
  void setup_test_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n";
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(UpdateNgpMesh, lazyAutoUpdate)
{
  setup_test_mesh();

  // Don't store persistent pointers/references if you want automatic updates
  // when acquiring an NgpMesh from BulkData
  stk::mesh::NgpMesh * ngpMesh = &get_bulk().get_ngp_mesh();

  get_bulk().modification_begin();
  get_bulk().modification_end();

#ifdef KOKKOS_ENABLE_CUDA
  EXPECT_FALSE(ngpMesh->is_up_to_date());
  ngpMesh = &get_bulk().get_ngp_mesh();  // Trigger update
  EXPECT_TRUE(ngpMesh->is_up_to_date());
#else
  EXPECT_TRUE(ngpMesh->is_up_to_date());
  ngpMesh = &get_bulk().get_ngp_mesh();
  EXPECT_TRUE(ngpMesh->is_up_to_date());
#endif
}

TEST_F(UpdateNgpMesh, manualUpdate)
{
  setup_test_mesh();

  // If storing a persistent reference, call the update_ngp_mesh() method
  // to ensure that it is synchronized with BulkData
  stk::mesh::NgpMesh & ngpMesh = get_bulk().get_ngp_mesh();

  get_bulk().modification_begin();
  get_bulk().modification_end();

#ifdef KOKKOS_ENABLE_CUDA
  EXPECT_FALSE(ngpMesh.is_up_to_date());
  get_bulk().update_ngp_mesh();
  EXPECT_TRUE(ngpMesh.is_up_to_date());
#else
  EXPECT_TRUE(ngpMesh.is_up_to_date());
  get_bulk().update_ngp_mesh();
  EXPECT_TRUE(ngpMesh.is_up_to_date());
#endif

}

}

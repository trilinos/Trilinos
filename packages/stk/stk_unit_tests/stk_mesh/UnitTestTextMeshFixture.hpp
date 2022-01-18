#ifndef STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_
#define STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>  // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/TextMeshFixture.hpp>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "mpi.h"

namespace
{
class TestTextMesh : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

class TestTextMeshAura : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMeshAura() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

class TestTextMesh2d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh2d() : TextMeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

class TestTextMesh1d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh1d() : TextMeshFixture(1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};
}



#endif /* STK_STK_UNIT_TESTS_STK_MESH_UNITTESTTEXTMESHFIXTURE_HPP_ */

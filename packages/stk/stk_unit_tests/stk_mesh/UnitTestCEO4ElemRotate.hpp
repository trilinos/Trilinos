#ifndef STK_STK_MESH_UNIT_TESTS_UNITTESTCEO4ELEMROTATE_HPP_
#define STK_STK_MESH_UNIT_TESTS_UNITTESTCEO4ELEMROTATE_HPP_

#include "stk_mesh/base/Types.hpp"

namespace stk { namespace unit_test_util { class BulkDataTester; } }
namespace stk { namespace mesh { class MetaData; } }

namespace CEOUtils
{

//////////////////////////////////// 4Elem4ProcRotate //////////////////////////////////////

void fillMeshfor4Elem4ProcRotateAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEOME_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

}

#endif

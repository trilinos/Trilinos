#ifndef STK_STK_MESH_UNIT_TESTS_UNITTESTCEO8ELEM_HPP_
#define STK_STK_MESH_UNIT_TESTS_UNITTESTCEO8ELEM_HPP_

#include "stk_mesh/base/Types.hpp"

namespace stk { namespace unit_test_util { class BulkDataTester; } }
namespace stk { namespace mesh { class MetaData; } }

namespace CEOUtils
{

//////////////////////////////////// 8Elem4ProcMoveTop //////////////////////////////////////

void fillMeshfor8Elem4ProcMoveTopAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh);

}

#endif

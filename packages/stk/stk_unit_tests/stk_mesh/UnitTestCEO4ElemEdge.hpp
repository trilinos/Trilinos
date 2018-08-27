#ifndef STK_STK_MESH_UNIT_TESTS_UNITTESTCEO4ELEMEDGE_HPP_
#define STK_STK_MESH_UNIT_TESTS_UNITTESTCEO4ELEMEDGE_HPP_

#include "stk_mesh/base/Types.hpp"

namespace stk { namespace unit_test_util { class BulkDataTester; } }
namespace stk { namespace mesh { class MetaData; } }

namespace CEOUtils
{

//////////////////////////////////// 4Elem4ProcEdge //////////////////////////////////////

void fillMeshfor4Elem4ProcEdgeAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data,
                                      stk::mesh::EntityKey &elem_key_chg_own);

void checkStatesAfterCEO_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh);

}

#endif

#ifndef STK_STK_MESH_UNIT_TESTS_UNITTESTCEO2ELEM_HPP_
#define STK_STK_MESH_UNIT_TESTS_UNITTESTCEO2ELEM_HPP_

#include "stk_unit_test_utils/BulkDataTester.hpp"           // for BulkDataTester
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include <vector>

namespace stk { namespace mesh { class MetaData; } }


namespace CEOUtils
{

//////////////////////////////////// 2Elem2ProcMove //////////////////////////////////////
void fillMeshfor2Elem2ProcMoveAndTest(stk::unit_test_util::BulkDataTester& bulk, stk::mesh::MetaData &meta, std::vector<stk::mesh::Entity>& elems);

void checkStatesAfterCEO_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk);

void checkStatesAfterCEOME_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk);

//////////////////////////////////// 2Elem2ProcFlip //////////////////////////////////////

void fillMeshfor2Elem2ProcFlipAndTest(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh);

void checkStatesAfterCEOME_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh);

//these tests are for turning regenerate_aura off in various places

void checkStatesAfterCEOME_2Elem2ProcMove_no_ghost(stk::unit_test_util::BulkDataTester &bulk);

void fillMeshfor2Elem2ProcFlipAndTest_no_ghost(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEOME_2Elem2ProcFlip_no_ghost(stk::unit_test_util::BulkDataTester& mesh);

}

#endif /* STK_STK_MESH_UNIT_TESTS_UNITTESTCEO2ELEM_HPP_ */

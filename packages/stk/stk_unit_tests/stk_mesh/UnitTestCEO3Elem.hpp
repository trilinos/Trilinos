#ifndef STK_STK_MESH_UNIT_TESTS_UNITTESTCEO3ELEM_HPP_
#define STK_STK_MESH_UNIT_TESTS_UNITTESTCEO3ELEM_HPP_

namespace stk { namespace unit_test_util { class BulkDataTester; } }
namespace stk { namespace mesh { class MetaData; } }

#include "stk_mesh/base/Types.hpp"

namespace CEOUtils
{

//////////////////////////////////// 3Elem2ProcMoveRight //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveRightAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector& elements);

void checkStatesAfterCEO_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 3Elem2ProcMoveLeft //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveLeftAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector &elements);

void checkStatesAfterCEO_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 3Elem4Proc1Edge3D //////////////////////////////////////

void fillMeshfor3Elem4Proc1Edge3DAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh);

}


#endif

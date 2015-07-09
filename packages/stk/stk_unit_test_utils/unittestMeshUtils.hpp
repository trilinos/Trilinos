#ifndef unittestMeshUtils_hpp
#define unittestMeshUtils_hpp


namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk
{
namespace unit_test_util
{


void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part);




} // namespace unit_test_util
} // namespace stk


#endif // unittestMeshUtils_hpp

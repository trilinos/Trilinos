#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace balance
{

stk::mesh::EntityIdProcMap make_mesh_consistent_with_parallel_mesh_rule1(stk::mesh::BulkData& bulkData);

}}

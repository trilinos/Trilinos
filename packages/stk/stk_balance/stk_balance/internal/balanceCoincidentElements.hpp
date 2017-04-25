namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace balance { class DecompositionChangeList; }}

namespace stk
{
namespace balance
{

void keep_coincident_elements_together(stk::mesh::BulkData &bulk, DecompositionChangeList &changeList);

}}

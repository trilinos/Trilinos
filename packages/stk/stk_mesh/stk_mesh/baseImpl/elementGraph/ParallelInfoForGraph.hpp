#ifndef STKMESHPARALLELINFOFORGRAPH_HPP
#define STKMESHPARALLELINFOFORGRAPH_HPP

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class ElemElemGraph; } }

namespace stk
{
namespace mesh
{
namespace impl
{

void update_parallel_graph_for_part_ordinals(ElemElemGraph& graph, const stk::mesh::BulkData& bulkData);

}
}
}

#endif

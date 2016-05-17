#ifndef STKMESHPARALLELINFOFORGRAPH_HPP
#define STKMESHPARALLELINFOFORGRAPH_HPP
#include <map>
#include <vector>
#include "../../base/Types.hpp"
#include "GraphTypes.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class ElemElemGraph; } }

namespace stk
{
namespace mesh
{
namespace impl
{

void populate_part_ordinals_for_remote_edges(const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo);
void populate_selected_value_for_remote_elements(const stk::mesh::BulkData& bulkData, ElemElemGraph& graph, stk::mesh::Selector selector, ParallelSelectedInfo &selInfo);

}
}
}

#endif

#ifndef STKMESHPARALLELINFOFORGRAPH_HPP
#define STKMESHPARALLELINFOFORGRAPH_HPP
#include <map>
#include <vector>
#include "../../base/Types.hpp"
#include "ElemElemGraphImpl.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class ElemElemGraph; } }

namespace stk
{
namespace mesh
{
namespace impl
{

typedef std::vector<stk::mesh::PartOrdinal> PartOrdinals;
typedef std::map<LocalId, PartOrdinals> ParallelPartInfo;

void populate_part_ordinals_for_remote_edges(const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo);
void update_parallel_graph_for_air_selector(const stk::mesh::BulkData& bulkData, ElemElemGraph& graph, stk::mesh::Selector air);
void update_parallel_graph_for_skin_selector(const stk::mesh::BulkData& bulkData, ElemElemGraph& graph, stk::mesh::Selector skinSelector);

}
}
}

#endif

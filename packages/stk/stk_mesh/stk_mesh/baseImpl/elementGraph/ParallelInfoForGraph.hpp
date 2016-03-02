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
typedef std::map<stk::mesh::GraphEdge, PartOrdinals, GraphEdgeLessByElem2> ParallelPartInfo;

void populate_part_ordinals_for_remote_edges(const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo);

}
}
}

#endif

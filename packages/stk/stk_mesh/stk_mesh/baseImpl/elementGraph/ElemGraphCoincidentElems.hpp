#ifndef STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP
#define STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/SparseGraph.hpp>

namespace stk { namespace mesh { class Graph; } }
namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }

namespace stk
{
namespace mesh
{
namespace impl
{

class IdMapper
{
public:
    virtual ~IdMapper() { }
    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId local) const = 0;
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const = 0;
};

bool is_coincident_connection(const stk::mesh::BulkData &bulkData,
                              stk::mesh::Entity localElem,
                              const stk::mesh::EntityVector& localElemSideNodes,
                              unsigned sideIndex,
                              stk::topology otherElemTopology,
                              const stk::mesh::EntityVector &otherElemSideNodes);

}}} // end namespaces stk mesh

#endif

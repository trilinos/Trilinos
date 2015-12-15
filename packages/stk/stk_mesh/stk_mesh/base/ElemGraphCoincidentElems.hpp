#ifndef STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP
#define STK_ELEM_ELEM_GRAPH_COINCIDENT_ELEMS_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class Graph; } }
namespace stk { namespace mesh { class ParallelInfoForGraphEdges; } }

namespace stk
{
namespace mesh
{
namespace impl
{

struct CoincidentElementDescription
{
    int numSides;
    stk::mesh::impl::LocalId elem1;
    stk::mesh::impl::LocalId elem2;
};

typedef std::map<stk::mesh::impl::LocalId, std::vector<stk::mesh::GraphEdge>> SparseGraph;

class IdMapper
{
public:
    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId local) const = 0;
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const = 0;
};

SparseGraph extract_coincident_sides(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies);
void append_extracted_coincident_sides(stk::mesh::Graph &graph,
                                       const std::vector<stk::topology> &topologies,
                                       const std::vector<impl::LocalId> &elemIds,
                                       stk::mesh::impl::SparseGraph &coincidentEdges);
void choose_face_id_for_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parallelInfoForGraphEdges,
                                            const stk::mesh::impl::SparseGraph &extractedCoincidentElements,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm);

}}} // end namespaces stk mesh

#endif

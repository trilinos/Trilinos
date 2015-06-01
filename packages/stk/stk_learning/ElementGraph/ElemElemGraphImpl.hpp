#ifndef STK_ELEM_ELEM_GRAPH_IMPL_HPP
#define STK_ELEM_ELEM_GRAPH_IMPL_HPP

#include <vector>
#include <map>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { class CommSparse; }

typedef int64_t LocalId;

//BeginDocExample3
struct parallel_info
{
    int m_other_proc;
    int m_other_side_ord;
    int m_permutation;
    bool m_in_part;
    parallel_info(int proc, int side_ord, int perm) :
        m_other_proc(proc), m_other_side_ord(side_ord), m_permutation(perm), m_in_part(true) {}
};
//EndDocExample3

typedef std::pair<LocalId,int> ElementSidePair;
typedef std::map<std::pair<LocalId,stk::mesh::EntityId>, parallel_info > ParallelGraphInfo;
typedef std::vector<std::vector<LocalId> > ElementGraph;
typedef std::vector<std::vector<int> > SidesForElementGraph;

namespace impl
{
void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies);

void fill_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides);

stk::mesh::EntityVector get_elements_to_communicate(const stk::mesh::BulkData& bulkData);

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate);

void add_possibly_connected_elements_to_graph_using_side_nodes(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, const stk::mesh::EntityVector& side_nodes, ParallelGraphInfo& parallel_graph_info,
        LocalId other_element, int other_side, int other_proc);

void fill_parallel_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info);
}

#endif

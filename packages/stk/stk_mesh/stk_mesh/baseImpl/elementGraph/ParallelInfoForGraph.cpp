#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include <stk_util/parallel/CommSparse.hpp>
#include "../EquivalentEntityBlocks.hpp"

namespace stk
{
namespace mesh
{
namespace impl
{
void pack_data_for_part_ordinals(stk::CommSparse &comm, ElemElemGraph& graph, const stk::mesh::BulkData& bulkData);
void pack_edge(stk::CommSparse &comm, ElemElemGraph& graph, const stk::mesh::BulkData& bulkData, const stk::mesh::GraphEdge& edge, int other_proc);
void unpack_and_update_part_ordinals(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, ElemElemGraph& graph);
stk::mesh::GraphEdge unpack_edge(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, ElemElemGraph& graph, int proc_id);

void update_parallel_graph_for_part_ordinals(ElemElemGraph& graph, const stk::mesh::BulkData& bulkData)
{
    stk::CommSparse comm(bulkData.parallel());
    pack_data_for_part_ordinals(comm, graph, bulkData);
    comm.allocate_buffers();
    pack_data_for_part_ordinals(comm, graph, bulkData);
    comm.communicate();
    unpack_and_update_part_ordinals(comm, bulkData, graph);
}

void pack_data_for_part_ordinals(stk::CommSparse &comm, ElemElemGraph& graph, const stk::mesh::BulkData& bulkData)
{
    const stk::mesh::impl::ParallelGraphInfo& parallel_info = graph.get_parallel_graph().get_parallel_graph_info();
    for(const auto& item : parallel_info)
    {
        const stk::mesh::GraphEdge &edge = item.first;
        const stk::mesh::impl::ParallelInfo &pinfo = item.second;
        stk::mesh::Entity local_element = graph.get_entity(edge.elem1);
        std::vector<stk::mesh::PartOrdinal> partOrdinals = stk::mesh::impl::get_element_block_part_ordinals(local_element, bulkData);

        bool didPartChangesOccur = true;

        if(didPartChangesOccur)
        {
            pack_edge(comm, graph, bulkData, edge, pinfo.get_proc_rank_of_neighbor());

            comm.send_buffer(pinfo.get_proc_rank_of_neighbor()).pack<size_t>(partOrdinals.size());
            for(stk::mesh::PartOrdinal partOrdinal : partOrdinals)
                comm.send_buffer(pinfo.get_proc_rank_of_neighbor()).pack<stk::mesh::PartOrdinal>(partOrdinal);
        }
    }
}

void pack_edge(stk::CommSparse &comm, ElemElemGraph& graph, const stk::mesh::BulkData& bulkData, const stk::mesh::GraphEdge& edge, int other_proc)
{
    stk::mesh::EntityId id1 = bulkData.identifier(graph.get_entity(edge.elem1));
    unsigned side1 = edge.side1;
    stk::mesh::EntityId id2 = -edge.elem2;
    unsigned side2 = edge.side2;
    comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(id1);
    comm.send_buffer(other_proc).pack<unsigned>(side1);
    comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(id2);
    comm.send_buffer(other_proc).pack<unsigned>(side2);
}

void unpack_and_update_part_ordinals(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, ElemElemGraph& graph)
{
    for(int i=0;i<bulkData.parallel_size();++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::GraphEdge edge = unpack_edge(comm, bulkData, graph, i);

            size_t num_ordinals = 0;
            comm.recv_buffer(i).unpack<size_t>(num_ordinals);
            std::vector<stk::mesh::PartOrdinal> partOrdinals(num_ordinals);
            for(stk::mesh::PartOrdinal &partOrdinal : partOrdinals)
                comm.recv_buffer(i).unpack<stk::mesh::PartOrdinal>(partOrdinal);

            stk::mesh::impl::ParallelGraphInfo& parallel_info = graph.get_parallel_graph().get_parallel_graph_info();
            auto iter = parallel_info.find(edge);
            stk::ThrowRequireWithSierraHelpMsg(iter!=parallel_info.end());
            iter->second.set_part_ordinals(partOrdinals);
        }
    }
}

stk::mesh::GraphEdge unpack_edge(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, ElemElemGraph& graph, int proc_id)
{
    stk::mesh::EntityId id1 = 0, id2 = 0;
    unsigned side1 = 0, side2 = 0;
    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(id1);
    comm.recv_buffer(proc_id).unpack<unsigned>(side1);
    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(id2);
    comm.recv_buffer(proc_id).unpack<unsigned>(side2);

    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, id2);
    stk::ThrowRequireWithSierraHelpMsg(bulkData.is_valid(element));

    stk::mesh::impl::LocalId localId2 = graph.get_local_element_id(element);
    stk::mesh::GraphEdge edge(localId2, side2, -id1, side1);
    return edge;
}

}
}
}

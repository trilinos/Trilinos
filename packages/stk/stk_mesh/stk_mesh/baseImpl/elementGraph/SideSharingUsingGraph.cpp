#include "SideSharingUsingGraph.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace stk { namespace mesh {

stk::mesh::EntityVector fill_shared_entities_that_need_fixing(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()), sides);

    stk::mesh::EntityVector sidesThatNeedFixing;
    for(stk::mesh::Entity side : sides)
        if(bulkData.state(side) == stk::mesh::Created)
        {
            unsigned num_nodes = bulkData.num_nodes(side);
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(side);

            std::vector<stk::mesh::EntityKey> nodeKeys(num_nodes);
            for(unsigned int i=0;i<num_nodes;++i)
                nodeKeys[i] = bulkData.entity_key(nodes[i]);

            std::vector<int> shared_procs;
            bulkData.shared_procs_intersection(nodeKeys, shared_procs);
            if(!shared_procs.empty())
                sidesThatNeedFixing.push_back(side);
        }
    return sidesThatNeedFixing;
}

stk::mesh::impl::ElementViaSidePair find_element_side_ord_for_side(const stk::mesh::BulkData& bulkData, unsigned num_elements, const stk::mesh::Entity* elements, stk::mesh::Entity side)
{
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    for(unsigned i=0;i<num_elements;++i)
    {
        if(bulkData.bucket(elements[i]).owned())
        {
            const stk::mesh::Entity* sides = bulkData.begin(elements[i], meta.side_rank());
            unsigned num_sides = bulkData.num_connectivity(elements[i], meta.side_rank());

            for(unsigned j=0;j<num_sides;++j)
            {
                if(sides[j]==side)
                {
                    const stk::mesh::ConnectivityOrdinal* ordinals = bulkData.begin_ordinals(elements[i], meta.side_rank());
                    return stk::mesh::impl::ElementViaSidePair{elements[i], static_cast<int>(ordinals[j])};
                }
            }
        }
    }
    return stk::mesh::impl::ElementViaSidePair{stk::mesh::Entity(),0};
}

stk::mesh::impl::ElementViaSidePair get_element_and_side_ordinal(const stk::mesh::BulkData& bulkData, stk::mesh::Entity side)
{
    unsigned num_elements = bulkData.num_elements(side);
    const stk::mesh::Entity* elements = bulkData.begin_elements(side);
    return find_element_side_ord_for_side(bulkData, num_elements, elements, side);
}

const unsigned * get_first_part_after_owned_and_shared(const stk::mesh::PartOrdinal sharedOrd, const unsigned *start, const unsigned *end)
{
    const unsigned skipUniversalAndOwned = 2;
    start += skipUniversalAndOwned;
    if(start < end && *start == sharedOrd)
        start++;
    return start;
}

void fill_part_ordinals_besides_owned_and_shared(const stk::mesh::Bucket &bucket, const stk::mesh::PartOrdinal sharedOrd, stk::mesh::OrdinalVector &partOrdinals)
{
    std::pair<const unsigned *, const unsigned *> partOrdinalRange = bucket.superset_part_ordinals();
    const unsigned *start = get_first_part_after_owned_and_shared(sharedOrd, partOrdinalRange.first, partOrdinalRange.second);
    partOrdinals.assign(start, partOrdinalRange.second);
}

void fill_sharing_data(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph &graph, const stk::mesh::EntityVector& sidesThatNeedFixing, std::vector<SideSharingData>& sideSharingDataThisProc, std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides)
{
    // Element 1, side 5: face 15
    // Element 2, side 3: face 23
    // Are these faces the same? Yes: delete face 23, then connect face 15 to element 2 with negative permutation

    const stk::mesh::PartOrdinal sharedOrd = bulkData.mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal();
    for(size_t i=0;i<sidesThatNeedFixing.size();++i)
    {
        stk::mesh::impl::ElementViaSidePair elementAndSide = get_element_and_side_ordinal(bulkData, sidesThatNeedFixing[i]);
        stk::mesh::impl::LocalId localElemId = graph.get_local_element_id(elementAndSide.element);
        for(const stk::mesh::GraphEdge& edge : graph.get_edges_for_element(localElemId))
        {
            if(edge.side1() == elementAndSide.side && edge.elem2() < 0)
            {
                const stk::mesh::impl::ParallelInfo &pInfo = graph.get_parallel_info_for_graph_edge(edge);
                const stk::mesh::Entity* nodes = bulkData.begin_nodes(sidesThatNeedFixing[i]);
                unsigned numNodes = bulkData.num_nodes(sidesThatNeedFixing[i]);

                SideSharingData localTemp({bulkData.identifier(elementAndSide.element), elementAndSide.side},
                                          sidesThatNeedFixing[i],
                                          pInfo.get_proc_rank_of_neighbor(),
                                          std::min(bulkData.parallel_rank(),pInfo.get_proc_rank_of_neighbor()),
                                          bulkData.identifier(sidesThatNeedFixing[i]));

                localTemp.sideNodes.resize(numNodes);
                for(unsigned j=0; j<numNodes; ++j) {
                    localTemp.sideNodes[j] = bulkData.identifier(nodes[j]);
                }

                fill_part_ordinals_besides_owned_and_shared(bulkData.bucket(sidesThatNeedFixing[i]), sharedOrd, localTemp.partOrdinals);

                sideSharingDataThisProc.push_back(localTemp);

                stk::mesh::EntityId localId = -edge.elem2();
                idAndSides.push_back({localId, edge.side2()});
            }
        }
    }
}

void pack_data(stk::CommSparse& comm, const std::vector<SideSharingData>& sideSharingDataThisProc, const std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides)
{
    for(size_t i=0;i<idAndSides.size();++i)
    {
        int other_proc = sideSharingDataThisProc[i].sharingProc;
        comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(idAndSides[i].id);
        comm.send_buffer(other_proc).pack<int>(idAndSides[i].side);
        comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(sideSharingDataThisProc[i].chosenSideId);
        comm.send_buffer(other_proc).pack<size_t>(sideSharingDataThisProc[i].sideNodes.size());
        for(stk::mesh::EntityId nodeId : sideSharingDataThisProc[i].sideNodes) {
            comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(nodeId);
        }
        comm.send_buffer(other_proc).pack<size_t>(sideSharingDataThisProc[i].partOrdinals.size());
        for(stk::mesh::PartOrdinal partOrd : sideSharingDataThisProc[i].partOrdinals) {
            comm.send_buffer(other_proc).pack<stk::mesh::PartOrdinal>(partOrd);
        }
    }
}

void allocate_and_send(stk::CommSparse& comm, const std::vector<SideSharingData>& sideSharingDataThisProc, const std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides)
{
    pack_data(comm, sideSharingDataThisProc, idAndSides);
    comm.allocate_buffers();
    pack_data(comm, sideSharingDataThisProc, idAndSides);
    comm.communicate();
}

template <typename T>
void unpack_vector(CommBuffer &buffer, std::vector<T> &vec)
{
    size_t size;
    buffer.unpack<size_t>(size);
    vec.resize(size);
    for(unsigned j=0; j<size; ++j) {
        buffer.unpack<T>(vec[j]);
    }
}

void unpack_data(stk::CommSparse& comm, int my_proc_id, int num_procs, std::vector<SideSharingData>& sideSharingDataThisProc)
{
    for(int i=0;i<num_procs;++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::impl::IdViaSidePair receivedIdSide;
            stk::mesh::EntityId chosenId;
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(receivedIdSide.id);
            comm.recv_buffer(i).unpack<int>(receivedIdSide.side);
            comm.recv_buffer(i).unpack<stk::mesh::EntityId>(chosenId);
            SideSharingData localTemp;
            localTemp.elementAndSide = receivedIdSide;
            localTemp.owningProc = std::min(my_proc_id, i);;
            localTemp.sharingProc = i;
            localTemp.chosenSideId = chosenId;
            unpack_vector(comm.recv_buffer(i), localTemp.sideNodes);
            unpack_vector(comm.recv_buffer(i), localTemp.partOrdinals);

            sideSharingDataThisProc.push_back(localTemp);
        }
    }
}

}
}

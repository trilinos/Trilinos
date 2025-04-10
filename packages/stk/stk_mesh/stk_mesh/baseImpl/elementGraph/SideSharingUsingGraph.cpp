#include "SideSharingUsingGraph.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

namespace stk { namespace mesh {

stk::mesh::EntityVector fill_shared_entities_that_need_fixing(const stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank)
{
    stk::mesh::EntityVector sides;
    const bool sortByGlobalId = false;
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    stk::mesh::Selector owned = meta.locally_owned_part();

    stk::mesh::EntityVector sidesThatNeedFixing;
    stk::mesh::EntityVector nodeVec;
    std::vector<int> shared_procs;

    stk::mesh::get_entities(bulkData, rank, owned, sides, sortByGlobalId);

    for(stk::mesh::Entity side : sides) {
        if(bulkData.state(side) == stk::mesh::Created)
        {
            unsigned num_nodes = bulkData.num_nodes(side);
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(side);

            nodeVec.resize(num_nodes);
            for(unsigned int i=0;i<num_nodes;++i)
                nodeVec[i] = nodes[i];

            bulkData.shared_procs_intersection(nodeVec, shared_procs);
            if(!shared_procs.empty()) {
                sidesThatNeedFixing.push_back(side);
            }
        }
    }

    return sidesThatNeedFixing;
}

stk::mesh::impl::ElementViaSidePair find_element_side_ord_for_side(const stk::mesh::BulkData& bulkData,
                                                                   unsigned num_elements,
                                                                   const stk::mesh::Entity* elements,
                                                                   stk::mesh::Entity side,
                                                                   std::vector<stk::mesh::ElemSideOrdinal>& scratchElemSideOrdinals)
{
    for(unsigned i=0;i<num_elements;++i)
    {
        stk::mesh::Entity elem = elements[i];

        if(bulkData.bucket(elem).owned())
        {
            stk::mesh::fill_element_side_entries(bulkData, elem, false, scratchElemSideOrdinals);

            for(const auto& entry : scratchElemSideOrdinals)
            {
                if(entry.side == side)
                {
                    return stk::mesh::impl::ElementViaSidePair{elem, static_cast<int>(entry.ordinal)};
                }
            }
        }
    }
    return stk::mesh::impl::ElementViaSidePair{stk::mesh::Entity(),0};
}

stk::mesh::impl::ElementViaSidePair get_element_and_side_ordinal(const stk::mesh::BulkData& bulkData, stk::mesh::Entity side,
                                                                 std::vector<stk::mesh::ElemSideOrdinal>& scratchElemSideOrdinals)
{
    unsigned num_elements = bulkData.num_elements(side);
    const stk::mesh::Entity* elements = bulkData.begin_elements(side);
    return find_element_side_ord_for_side(bulkData, num_elements, elements, side, scratchElemSideOrdinals);
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

    std::vector<stk::mesh::ElemSideOrdinal> scratchElemSideOrdinals;
    const stk::mesh::PartOrdinal sharedOrd = bulkData.mesh_meta_data().globally_shared_part().mesh_meta_data_ordinal();
    for(size_t i=0;i<sidesThatNeedFixing.size();++i)
    {
        stk::mesh::impl::ElementViaSidePair elementAndSide = get_element_and_side_ordinal(bulkData, sidesThatNeedFixing[i], scratchElemSideOrdinals);
        if (elementAndSide.element.local_offset() == 0) {
            //solo side !!!
            continue;
        }
        stk::mesh::impl::LocalId localElemId = graph.get_local_element_id(elementAndSide.element);

        std::vector<int> sharedCoincProcs;
        int minSharedCoincProc = bulkData.parallel_rank();

        sharedCoincProcs.push_back(bulkData.parallel_rank());
        for(const stk::mesh::GraphEdge& edge : graph.get_edges_for_element(localElemId))
        {
            if(edge.side1() == elementAndSide.side && edge.elem2() < 0)
            {
                const stk::mesh::impl::ParallelInfo &pInfo = graph.get_parallel_info_for_graph_edge(edge);
                sharedCoincProcs.push_back( pInfo.get_proc_rank_of_neighbor() );
                minSharedCoincProc = std::min(minSharedCoincProc, pInfo.get_proc_rank_of_neighbor());

                stk::mesh::EntityId globalId = -edge.elem2();
                idAndSides.push_back({globalId, edge.side2()});
            }
        }


        SideSharingData localTemp( {bulkData.identifier(elementAndSide.element), elementAndSide.side},
                                  sidesThatNeedFixing[i],
                                  -1,
                                  sharedCoincProcs,
                                  minSharedCoincProc,
                                  bulkData.identifier(sidesThatNeedFixing[i]));

        fill_part_ordinals_besides_owned_and_shared(bulkData.bucket(sidesThatNeedFixing[i]), sharedOrd, localTemp.partOrdinals);

        for(int sharingProc : sharedCoincProcs)
        {
            if(sharingProc != bulkData.parallel_rank())
            {
                localTemp.sharingProc = sharingProc;
                sideSharingDataThisProc.push_back(localTemp);
            }
        }
    }
}

template <typename T>
void pack_vector(CommBuffer &buffer, const std::vector<T> &vec)
{
    buffer.pack<size_t>(vec.size());
    for(const T& entry : vec)
        buffer.pack<T>(entry);
}

void pack_data(stk::CommSparse& comm, const std::vector<SideSharingData>& sideSharingDataThisProc, const std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides)
{
    for(size_t i=0;i<idAndSides.size();++i)
    {
        int other_proc = sideSharingDataThisProc[i].sharingProc;
        CommBuffer &buffer = comm.send_buffer(other_proc);

        buffer.pack<stk::mesh::EntityId>(idAndSides[i].id);
        buffer.pack<int>(idAndSides[i].side);
        buffer.pack<stk::mesh::EntityId>(sideSharingDataThisProc[i].chosenSideId);
        pack_vector(buffer, sideSharingDataThisProc[i].allSharingProcs);
        pack_vector(buffer, sideSharingDataThisProc[i].partOrdinals);
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
            localTemp.owningProc = std::min(my_proc_id, i);
            localTemp.sharingProc = i;
            localTemp.chosenSideId = chosenId;
            unpack_vector(comm.recv_buffer(i), localTemp.allSharingProcs);
            unpack_vector(comm.recv_buffer(i), localTemp.partOrdinals);

            sideSharingDataThisProc.push_back(localTemp);
        }
    }
}

}
}

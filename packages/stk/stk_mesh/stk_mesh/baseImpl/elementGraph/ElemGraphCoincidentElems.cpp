#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "ElemElemGraph.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

void choose_consistent_face_ids_on_procs_that_own_coincident_elements(const stk::mesh::Graph& graph,
                                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                 const stk::mesh::impl::SparseGraph& extractedCoincidentElements,
                                                 const IdMapper &idMapper,
                                                 std::vector<stk::mesh::impl::ElementSidePair>& elemIdsToSendChangesFor);

void update_chosen_ids_on_other_procs_for_edges_with_coincident_elements(const std::vector<stk::mesh::impl::ElementSidePair> &elemIdsToSendChangesFor,
                                        const stk::mesh::Graph &graph,
                                        const stk::mesh::impl::SparseGraph &coincidentEdges,
                                        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                        const IdMapper& idMapper,
                                        MPI_Comm comm);



// AefB.e: A and B are ANCE, e is MCE, f is SCE ( e is MCE because it has smallest id of stacked elements e,f )
// master coincident element        MCE
// slave coincident element         SCE
// adjacent non-coincident element  ANCE

void make_chosen_ids_in_parinfo_consistent_for_edges_with_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parInfosForEdges,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm)
{
    std::vector<stk::mesh::impl::ElementSidePair> elemIdsToSendChangesFor;
    choose_consistent_face_ids_on_procs_that_own_coincident_elements(graph, parInfosForEdges, coincidentEdges, idMapper, elemIdsToSendChangesFor);
    update_chosen_ids_on_other_procs_for_edges_with_coincident_elements(elemIdsToSendChangesFor, graph, coincidentEdges, parInfosForEdges, idMapper, comm);
}


int count_shared_sides(const stk::mesh::Graph &graph, stk::mesh::impl::LocalId elem1, stk::mesh::impl::LocalId elem2)
{
    int numSharedSides = 0;
    for(size_t i=0; i < graph.get_num_edges_for_element(elem1); i++)
    {
        const stk::mesh::GraphEdge &graphEdge = graph.get_edge_for_element(elem1, i);
        if(graphEdge.elem2 == elem2)
            numSharedSides++;
    }
    return numSharedSides;
}

bool are_elements_coincident(const stk::mesh::Graph &graph,
                             const CoincidentElementDescription &elemDesc)
{
    int numSharedSides = count_shared_sides(graph, elemDesc.elem1, elemDesc.elem2);
    return (numSharedSides == elemDesc.numSides);
}

void fill_coincident_edges(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies, size_t elemId, std::vector<stk::mesh::GraphEdge> &edgesToDelete)
{
    for(size_t edgeIndex = 0; edgeIndex < graph.get_num_edges_for_element(elemId); edgeIndex++)
    {
        const stk::mesh::GraphEdge &graphEdge = graph.get_edge_for_element(elemId, edgeIndex);
        if(stk::mesh::impl::are_elements_coincident(graph, {static_cast<int>(topologies[elemId].num_sides()), graphEdge.elem1, graphEdge.elem2}))
            edgesToDelete.push_back(graphEdge);
    }
}

void delete_edges(stk::mesh::Graph& graph, const std::vector<stk::mesh::GraphEdge>& edgesToDelete)
{
    for(const stk::mesh::GraphEdge& edgeToDelete : edgesToDelete)
        graph.delete_edge(edgeToDelete);
}

void add_edges(const std::vector<stk::mesh::GraphEdge>& edgesToDelete, SparseGraph& extractedCoincidentSides)
{
    for(const stk::mesh::GraphEdge& edgeToDelete : edgesToDelete)
        extractedCoincidentSides[edgeToDelete.elem1].push_back(edgeToDelete);
}

void extract_coincident_sides_for_element(stk::mesh::Graph& graph,
                                          const std::vector<stk::topology>& topologies,
                                          size_t elemId,
                                          SparseGraph& extractedCoincidentSides)
{
    std::vector<stk::mesh::GraphEdge> coincidentEdges;
    fill_coincident_edges(graph, topologies, elemId, coincidentEdges);
    add_edges(coincidentEdges, extractedCoincidentSides);
    delete_edges(graph, coincidentEdges);
}

SparseGraph extract_coincident_sides(stk::mesh::Graph &graph, const std::vector<stk::topology> &topologies)
{
    SparseGraph extractedCoincidentSides;
    for(size_t elemId = 0; elemId < graph.get_num_elements_in_graph(); elemId++)
        extract_coincident_sides_for_element(graph, topologies, elemId, extractedCoincidentSides);
    return extractedCoincidentSides;
}

void append_extracted_coincident_sides(stk::mesh::Graph &graph,
                                       const std::vector<stk::topology> &topologies,
                                       const std::vector<impl::LocalId> &elemIds,
                                       stk::mesh::impl::SparseGraph &coincidentEdges)
{
    for(impl::LocalId elemId : elemIds)
        extract_coincident_sides_for_element(graph, topologies, elemId, coincidentEdges);
}

void pack_graph_edge_and_chosen_side_id_to_proc(stk::CommSparse& commSparse, int otherProc, const stk::mesh::GraphEdge& graphEdge, stk::mesh::EntityId chosenSideId)
{
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(graphEdge.elem1);
    commSparse.send_buffer(otherProc).pack<int>(graphEdge.side1);
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(graphEdge.elem2);
    commSparse.send_buffer(otherProc).pack<int>(graphEdge.side2);
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(chosenSideId);
}

typedef std::map<stk::mesh::GraphEdge, stk::mesh::impl::parallel_info *, GraphEdgeLessByElem2> GraphEdgeToParInfoMap;

void pack_graph_edge_and_chosen_side_id(stk::CommSparse& commSparse,
                                        const stk::mesh::GraphEdge& graphEdge,
                                        const parallel_info* parInfo,
                                        const IdMapper& idMapper)
{
    stk::mesh::GraphEdge remoteEdge(-graphEdge.elem2, graphEdge.side2, idMapper.localToGlobal(graphEdge.elem1), graphEdge.side1);
    pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfo->m_other_proc, remoteEdge, parInfo->m_chosen_side_id);
}


bool is_graph_edge_for_elem_side(const stk::mesh::GraphEdge& graphEdge, stk::mesh::impl::LocalId elemId, int side)
{
    return graphEdge.elem1 == elemId && graphEdge.side1 == side;
}

void pack_changed_edges_and_chosen_side_ids(stk::CommSparse& commSparse,
                                            const std::vector<stk::mesh::impl::ElementSidePair> &elemIdsToSendChangesFor,
                                            const stk::mesh::Graph &graph,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                            const IdMapper& idMapper)
{
    //for coincident side e,s that owns the choice of side id
    for(const stk::mesh::impl::ElementSidePair &elemSidePair : elemIdsToSendChangesFor)
    {
        for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(elemSidePair.first))
        {
            //is adjacent graph edge for e,s
            if(is_graph_edge_for_elem_side(graphEdge, elemSidePair.first, elemSidePair.second))
            {
                //if adjacent edge is remote
                if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
                {
                    //send edge e,s to adjacent elem to other proc with updated chosen side id
                    const stk::mesh::impl::parallel_info &parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                    pack_graph_edge_and_chosen_side_id(commSparse, graphEdge, &parInfo, idMapper);

                    // loop over elements coincident for e,s
                    auto iter = coincidentEdges.find(elemSidePair.first);
                    for(const stk::mesh::GraphEdge &coincidentEdge : iter->second)
                    {
                        // if coincident elem f is on another proc
                        if(is_graph_edge_for_elem_side(coincidentEdge, elemSidePair.first, elemSidePair.second) && !stk::mesh::impl::is_local_element(coincidentEdge.elem2))
                        {
                            //send edge from f,s to adjacent elem to other proc
                            const stk::mesh::impl::parallel_info &parInfoCoincident = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentEdge);
                            bool isAdjacentElementAndCoincidentElementBothOnSameProc = parInfoCoincident.m_other_proc == parInfo.m_other_proc;
                            if(!isAdjacentElementAndCoincidentElementBothOnSameProc)
                            {
                                stk::mesh::GraphEdge remoteEdge(-graphEdge.elem2, graphEdge.side2, -coincidentEdge.elem2, coincidentEdge.side2);
                                pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfo.m_other_proc, remoteEdge, parInfo.m_chosen_side_id);
                                stk::mesh::GraphEdge remoteEdge2(-coincidentEdge.elem2, coincidentEdge.side2, -graphEdge.elem2, graphEdge.side2);
                                pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoCoincident.m_other_proc, remoteEdge2, parInfoCoincident.m_chosen_side_id);
                            }
                        }
                    }
                }
                //adjacent elem is local
                else
                {
                    // loop over elements coincident for e,s
                    auto iter = coincidentEdges.find(elemSidePair.first);
                    for(const stk::mesh::GraphEdge &coincidentEdge : iter->second)
                    {
                        // if coincident elem f is on another proc
                        if(is_graph_edge_for_elem_side(coincidentEdge, elemSidePair.first, elemSidePair.second) && !stk::mesh::impl::is_local_element(coincidentEdge.elem2))
                        {
                            const stk::mesh::impl::parallel_info &parInfoCoincident = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentEdge);
                            stk::mesh::GraphEdge remoteEdge2(-coincidentEdge.elem2, coincidentEdge.side2, idMapper.localToGlobal(graphEdge.elem2), graphEdge.side2);
                            pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoCoincident.m_other_proc, remoteEdge2, parInfoCoincident.m_chosen_side_id);
                        }
                    }
                }
            }
        }

        //send coincident edge e,0->f,0 to proc that has f
        auto iter = coincidentEdges.find(elemSidePair.first);
        for(const stk::mesh::GraphEdge &coincidentEdge : iter->second)
        {
            if(is_graph_edge_for_elem_side(coincidentEdge, elemSidePair.first, elemSidePair.second) && !stk::mesh::impl::is_local_element(coincidentEdge.elem2))
            {
                const stk::mesh::impl::parallel_info &parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentEdge);
                pack_graph_edge_and_chosen_side_id(commSparse, coincidentEdge, &parInfo, idMapper);
            }
        }
    }
}

void pack_changed_edges_and_chosen_side_ids2(stk::CommSparse& commSparse,
                                            const std::vector<stk::mesh::impl::ElementSidePair> &elemIdsToSendChangesFor,
                                            const stk::mesh::Graph &graph,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                            const IdMapper& idMapper)
{
    //for coincident side e,s that owns the choice of side id
    for(const stk::mesh::impl::ElementSidePair &elemSidePair : elemIdsToSendChangesFor)
    {
        // loop over elements coincident for e,s
        auto iter = coincidentEdges.find(elemSidePair.first);
        for(const stk::mesh::GraphEdge &coincidentEdge : iter->second)
        {
            if(is_graph_edge_for_elem_side(coincidentEdge, elemSidePair.first, elemSidePair.second))
            {
                // if coincident elem f is on another proc
                if(!stk::mesh::impl::is_local_element(coincidentEdge.elem2))
                {
                    const stk::mesh::impl::parallel_info &parInfoCoincident = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentEdge);

                    //send coincident edge e,0->f,0 to proc that has f
                    pack_graph_edge_and_chosen_side_id(commSparse, coincidentEdge, &parInfoCoincident, idMapper);

                    for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(elemSidePair.first))
                    {
                        //is adjacent graph edge for e,s
                        if(is_graph_edge_for_elem_side(graphEdge, elemSidePair.first, elemSidePair.second))
                        {
                            //if adjacent edge is remote
                            if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
                            {
                                //send edge e,s->a,? to other proc with updated chosen side id
                                const stk::mesh::impl::parallel_info &parInfoAdjacent = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                                pack_graph_edge_and_chosen_side_id(commSparse, graphEdge, &parInfoAdjacent, idMapper);

                                //send edge from f,s to adjacent elem to other proc
                                bool isAdjacentElementAndCoincidentElementBothOnSameProc = parInfoCoincident.m_other_proc == parInfoAdjacent.m_other_proc;
                                if(!isAdjacentElementAndCoincidentElementBothOnSameProc)
                                {
                                    stk::mesh::GraphEdge remoteEdge(-graphEdge.elem2, graphEdge.side2, -coincidentEdge.elem2, coincidentEdge.side2);
                                    pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoAdjacent.m_other_proc, remoteEdge, parInfoAdjacent.m_chosen_side_id);
                                    stk::mesh::GraphEdge remoteEdge2(-coincidentEdge.elem2, coincidentEdge.side2, -graphEdge.elem2, graphEdge.side2);
                                    pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoCoincident.m_other_proc, remoteEdge2, parInfoCoincident.m_chosen_side_id);
                                }
                            }
                            //adjacent elem is local
                            else
                            {
                                const stk::mesh::impl::parallel_info &parInfoCoincident = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentEdge);
                                stk::mesh::GraphEdge remoteEdge2(-coincidentEdge.elem2, coincidentEdge.side2, idMapper.localToGlobal(graphEdge.elem2), graphEdge.side2);
                                pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoCoincident.m_other_proc, remoteEdge2, parInfoCoincident.m_chosen_side_id);
                            }
                        }
                    }
                }
                // coincident elem f is local
                else
                {
                    for(const stk::mesh::GraphEdge &graphEdge : graph.get_edges_for_element(elemSidePair.first))
                    {
                        //is adjacent graph edge for e,s
                        if(is_graph_edge_for_elem_side(graphEdge, elemSidePair.first, elemSidePair.second) && !stk::mesh::impl::is_local_element(graphEdge.elem2))
                        {
                            //send edge f,s->a,? to other proc with updated chosen side id
                            stk::mesh::GraphEdge adjToCoincidentEdge(-graphEdge.elem2, graphEdge.side2, idMapper.localToGlobal(coincidentEdge.elem2), coincidentEdge.side2);
                            const stk::mesh::impl::parallel_info &parInfoAdjacent = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(graphEdge);
                            pack_graph_edge_and_chosen_side_id_to_proc(commSparse, parInfoAdjacent.m_other_proc, adjToCoincidentEdge, parInfoAdjacent.m_chosen_side_id);
                        }
                    }
                }
            }
        }
    }
}

void unpack_local_elem_side(stk::CommSparse& commSparse,
                            int procId,
                            const IdMapper& idMapper,
                            stk::mesh::GraphEdge& recvGraphEdge)
{
    stk::mesh::EntityId globalId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(globalId);
    recvGraphEdge.elem1 = idMapper.globalToLocal(globalId);
    commSparse.recv_buffer(procId).unpack<int>(recvGraphEdge.side1);
}

void unpack_remote_elem_side(stk::CommSparse& commSparse, int procId, stk::mesh::GraphEdge& recvGraphEdge)
{
    stk::mesh::EntityId globalId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(globalId);
    recvGraphEdge.elem2 = -globalId;
    commSparse.recv_buffer(procId).unpack<int>(recvGraphEdge.side2);
}

stk::mesh::GraphEdge unpack_graph_edge(stk::CommSparse& commSparse, int procId, const IdMapper& idMapper)
{
    stk::mesh::GraphEdge recvGraphEdge;
    unpack_local_elem_side(commSparse, procId, idMapper, recvGraphEdge);
    unpack_remote_elem_side(commSparse, procId, recvGraphEdge);
    return recvGraphEdge;
}

void update_chosen_side_id_for_graph_edge(stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                          const stk::mesh::GraphEdge& recvGraphEdge,
                                          stk::mesh::EntityId chosenSideId)
{
    stk::mesh::impl::parallel_info& parInfoToChange = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(recvGraphEdge);
    parInfoToChange.m_chosen_side_id = chosenSideId;
}

stk::mesh::EntityId unpack_chosen_side_id(stk::CommSparse &commSparse, int procId)
{
    stk::mesh::EntityId chosenSideId;
    commSparse.recv_buffer(procId).unpack<stk::mesh::EntityId>(chosenSideId);
    return chosenSideId;
}

void unpack_changed_edges_and_chosen_side_ids(stk::CommSparse &commSparse,
                                              int procId,
                                              const IdMapper& idMapper,
                                              stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    stk::mesh::GraphEdge recvGraphEdge = unpack_graph_edge(commSparse, procId, idMapper);
    stk::mesh::EntityId chosenSideId = unpack_chosen_side_id(commSparse, procId);
    update_chosen_side_id_for_graph_edge(parallelInfoForGraphEdges, recvGraphEdge, chosenSideId);
}

void update_chosen_ids_on_other_procs_for_edges_with_coincident_elements(const std::vector<stk::mesh::impl::ElementSidePair> &elemIdsToSendChangesFor,
                                        const stk::mesh::Graph &graph,
                                        const stk::mesh::impl::SparseGraph &coincidentEdges,
                                        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                        const IdMapper& idMapper,
                                        MPI_Comm comm)
{
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&commSparse, &elemIdsToSendChangesFor, &graph, &coincidentEdges, &parallelInfoForGraphEdges, &idMapper]()
        {
            pack_changed_edges_and_chosen_side_ids2(commSparse, elemIdsToSendChangesFor, graph, coincidentEdges, parallelInfoForGraphEdges, idMapper);
        });
    stk::unpack_communications(commSparse,
        [&commSparse, &idMapper, &parallelInfoForGraphEdges](int procId)
        {
            unpack_changed_edges_and_chosen_side_ids(commSparse, procId, idMapper, parallelInfoForGraphEdges);
        });
}

void add_graph_edge_and_par_info(const stk::mesh::GraphEdge& parallelGraphEdge,
                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                 GraphEdgeToParInfoMap& parInfos)
{
    stk::mesh::impl::parallel_info& parInfoThisElem =
            parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(parallelGraphEdge);
    parInfos.insert(std::make_pair(parallelGraphEdge, &parInfoThisElem));
}

bool are_graph_edges_for_same_side(const stk::mesh::GraphEdge& parallelGraphEdge, const stk::mesh::GraphEdge& coincidentGraphEdge)
{
    return is_graph_edge_for_elem_side(parallelGraphEdge, coincidentGraphEdge.elem1, coincidentGraphEdge.side1);
}

void add_coincident_graph_edges_and_par_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                              const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                                              stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                              GraphEdgeToParInfoMap& parInfos)
{
    for(const stk::mesh::GraphEdge& coincidentGraphEdge : coincidentEdgesForElem)
    {
        if(are_graph_edges_for_same_side(parallelGraphEdge, coincidentGraphEdge))
        {
            if (impl::is_local_element(coincidentGraphEdge.elem2))
            {
                stk::mesh::GraphEdge otherElemGraphEdge(coincidentGraphEdge.elem2,
                                                        coincidentGraphEdge.side2,
                                                        parallelGraphEdge.elem2,
                                                        parallelGraphEdge.side2);
                // add (SCE, ANCE) edge to parInfo
                add_graph_edge_and_par_info(otherElemGraphEdge, parallelInfoForGraphEdges, parInfos);
            }
            else
            {
                // add (MCE, SCE) edge to parInfo
                add_graph_edge_and_par_info(coincidentGraphEdge, parallelInfoForGraphEdges, parInfos);
            }
        }
    }
}

GraphEdgeToParInfoMap get_coincident_parallel_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                                    const std::vector<stk::mesh::GraphEdge> &coincidentEdgesForElem,
                                                    stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    GraphEdgeToParInfoMap parInfos;
    // add (MCE, ANCE) edge to parInfo
    add_graph_edge_and_par_info(parallelGraphEdge, parallelInfoForGraphEdges, parInfos);
    // add (MCE, SCE) or (SCE, ANCE) edges to parInfo
    add_coincident_graph_edges_and_par_infos(parallelGraphEdge, coincidentEdgesForElem, parallelInfoForGraphEdges, parInfos);
    return parInfos;
}

stk::mesh::EntityId get_min_chosen_side_id(const GraphEdgeToParInfoMap &coincidentParInfos)
{

    stk::mesh::EntityId chosenSideId = std::numeric_limits<stk::mesh::EntityId>::max();
    for(const GraphEdgeToParInfoMap::value_type &graphEdgeAndParInfo : coincidentParInfos)
        chosenSideId = std::min(chosenSideId, graphEdgeAndParInfo.second->m_chosen_side_id);
    return chosenSideId;
}

void update_chosen_side_id_for_coincident_graph_edges(const GraphEdgeToParInfoMap &coincidentParInfos,
                                                      stk::mesh::EntityId chosenSideId)
{
    for(const GraphEdgeToParInfoMap::value_type &graphEdgeAndParInfo : coincidentParInfos)
        graphEdgeAndParInfo.second->m_chosen_side_id = chosenSideId;
}

void match_ids_for_edges_related_to_this_side(const stk::mesh::GraphEdge& MCEtoANCEedge,
                                                  stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                  const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem)
{
    // For a particular side ordinal, get all related edges
    // add (MCE, ANCE) edge to parInfo and add either (MCE, SCE(1,2,3...)) or (SCE(1,2,3...), ANCE) edges to parInfo
    GraphEdgeToParInfoMap coincidentParInfos = get_coincident_parallel_infos(MCEtoANCEedge,
                                                                             coincidentEdgesForElem,
                                                                             parallelInfoForGraphEdges);
    stk::mesh::EntityId chosenSideId = get_min_chosen_side_id(coincidentParInfos);

    update_chosen_side_id_for_coincident_graph_edges(coincidentParInfos, chosenSideId);
}


void match_chosen_ids_for_element_edges_involving_coincident_elements(const stk::mesh::Graph& graph,
        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
        const stk::mesh::impl::LocalId MCE,
        std::vector<stk::mesh::impl::ElementSidePair>& elemIdsToSendChangesFor)
{
    for(const stk::mesh::GraphEdge& MCEtoANCEedge : graph.get_edges_for_element(MCE))
    {
        if(!stk::mesh::impl::is_local_element(MCEtoANCEedge.elem2))
        {
            elemIdsToSendChangesFor.push_back(stk::mesh::impl::ElementSidePair(MCEtoANCEedge.elem1, MCEtoANCEedge.side1));
            match_ids_for_edges_related_to_this_side(MCEtoANCEedge,
                                                     parallelInfoForGraphEdges,
                                                     coincidentEdgesForElem);
        }
    }
}

void match_chosen_ids_for_ANCE_to_SCE_edges(const stk::mesh::Graph& graph,
        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
        std::vector<stk::mesh::impl::ElementSidePair>& elemIdsToSendChangesFor)
{
    for(const stk::mesh::GraphEdge& coincidentGraphEdge : coincidentEdgesForElem)
    {
        if(!stk::mesh::impl::is_local_element(coincidentGraphEdge.elem2))
        {
            stk::mesh::impl::parallel_info& parInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentGraphEdge);
            stk::mesh::EntityId chosenSideId = parInfo.m_chosen_side_id;
            for(const stk::mesh::GraphEdge& graphEdge : graph.get_edges_for_element(coincidentGraphEdge.elem1))
            {
                if(are_graph_edges_for_same_side(coincidentGraphEdge, graphEdge))
                {
                    if(stk::mesh::impl::is_local_element(graphEdge.elem2))
                    {
                        for(const stk::mesh::GraphEdge& coincidentGraphEdge2 : coincidentEdgesForElem)
                        {
                            if(are_graph_edges_for_same_side(coincidentGraphEdge, coincidentGraphEdge2))
                            {
                                stk::mesh::GraphEdge coincidentTransposedEdge(graphEdge.elem2, graphEdge.side2, coincidentGraphEdge2.elem2, coincidentGraphEdge2.side2);
                                stk::mesh::impl::parallel_info& parInfoThisElem = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(coincidentTransposedEdge);
                                // set chosen id for parallel edge (ANCE, SCE(1,2,...))
                                parInfoThisElem.m_chosen_side_id = chosenSideId;
                            }
                        }
                    }
                }
            }
            // communicate (e, side1) to other procs
            elemIdsToSendChangesFor.push_back(stk::mesh::impl::ElementSidePair(coincidentGraphEdge.elem1, coincidentGraphEdge.side1));
        }
    }
}


void match_chosen_ids_for_edges_this_proc(const stk::mesh::Graph& graph,
                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                            const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                            const stk::mesh::impl::LocalId MCE,
                            std::vector<stk::mesh::impl::ElementSidePair>& elemIdsToSendChangesFor)
{
    // For a particular side ordinal, get all related edges
    // set chosen ids for (MCE, ANCE) edge and (MCE, SCE(1,2,3...)) or (SCE(1,2,3...), ANCE) edges
    match_chosen_ids_for_element_edges_involving_coincident_elements(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, MCE, elemIdsToSendChangesFor);
    // set chosen id for parallel edge (ANCE, SCE(1,2,...))
    match_chosen_ids_for_ANCE_to_SCE_edges(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, elemIdsToSendChangesFor);
}

bool do_i_own_coincident_element(const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem, const stk::mesh::impl::LocalId elemId, const IdMapper &idMapper)
{
    for(const stk::mesh::GraphEdge& graphEdge : coincidentEdgesForElem)
        if(idMapper.localToGlobal(elemId) > idMapper.localToGlobal(graphEdge.elem2))
            return false;
    return true;
}

void choose_consistent_face_ids_on_procs_that_own_coincident_elements(const stk::mesh::Graph& graph,
                                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                 const stk::mesh::impl::SparseGraph& extractedCoincidentElements,
                                                 const IdMapper &idMapper,
                                                 std::vector<stk::mesh::impl::ElementSidePair>& elemIdsToSendChangesFor)
{
    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : extractedCoincidentElements)
    {
        const stk::mesh::impl::LocalId possibleMCE = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        if(do_i_own_coincident_element(coincidentEdgesForElem, possibleMCE, idMapper))
            match_chosen_ids_for_edges_this_proc(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, possibleMCE, elemIdsToSendChangesFor);
    }
}

}}} // end namespaces stk mesh impl


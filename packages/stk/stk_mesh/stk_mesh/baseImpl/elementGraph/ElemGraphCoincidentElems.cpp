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

bool are_side_nodes_degenerate(const stk::mesh::EntityVector &sideNodes)
{
    std::set<stk::mesh::Entity> sortedNodes;

    for(const stk::mesh::Entity& node : sideNodes)
        sortedNodes.insert(node);

    return sortedNodes.size() != sideNodes.size();
}

struct TopologyChecker
{
    bool are_both_shells() const
    {
        return localTopology.is_shell() && remoteTopology.is_shell();
    }

    bool are_both_not_shells() const
    {
        return !localTopology.is_shell() && !remoteTopology.is_shell();
    }

    stk::topology localTopology;
    stk::topology remoteTopology;
};

bool is_side_node_permutation_positive(const stk::mesh::BulkData &bulkData, stk::mesh::Entity localElem, unsigned sideIndex, const stk::mesh::EntityVector &otherElemSideNodes)
{
    stk::mesh::EntityVector thisElemSideNodes = impl::get_element_side_nodes_from_topology(bulkData, localElem, sideIndex);
    stk::topology localTopology = bulkData.bucket(localElem).topology();
    stk::topology localSideTopology = localTopology.side_topology(sideIndex);
    std::pair<bool, unsigned> result = localSideTopology.equivalent(otherElemSideNodes, thisElemSideNodes);
    return result.first && localSideTopology.is_positive_polarity(result.second);
}

bool is_nondegenerate_coincident_connection(const stk::mesh::BulkData &bulkData,
                                            stk::mesh::Entity localElem,
                                            unsigned sideIndex,
                                            stk::topology otherElemTopology,
                                            const stk::mesh::EntityVector &otherElemSideNodes)
{
    stk::topology localTopology = bulkData.bucket(localElem).topology();
    TopologyChecker topologyChecker {localTopology, otherElemTopology};
    if(topologyChecker.are_both_shells())
        return true;
    if(topologyChecker.are_both_not_shells())
        return is_side_node_permutation_positive(bulkData, localElem, sideIndex, otherElemSideNodes);
    return false;
}

bool is_coincident_connection(const stk::mesh::BulkData &bulkData,
                              stk::mesh::Entity localElem,
                              unsigned sideIndex,
                              stk::topology otherElemTopology,
                              const stk::mesh::EntityVector &otherElemSideNodes)
{
    if(are_side_nodes_degenerate(otherElemSideNodes))
        return false;
    return is_nondegenerate_coincident_connection(bulkData, localElem, sideIndex, otherElemTopology, otherElemSideNodes);
}

typedef std::map<stk::mesh::impl::ElementSidePair, std::vector<stk::mesh::GraphEdge>> ElemSideAndEdges;

// AefB.e: A and B are ANCE, e is MCE, f is SCE ( e is MCE because it has smallest id of stacked elements e,f )
// master coincident element        MCE
// slave coincident element         SCE
// adjacent non-coincident element  ANCE

void pack_graph_edge_and_chosen_side_id_to_proc(stk::CommSparse& commSparse, int otherProc, const stk::mesh::GraphEdge& graphEdge, stk::mesh::EntityId chosenSideId)
{
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(graphEdge.elem1);
    commSparse.send_buffer(otherProc).pack<int>(graphEdge.side1);
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(graphEdge.elem2);
    commSparse.send_buffer(otherProc).pack<int>(graphEdge.side2);
    commSparse.send_buffer(otherProc).pack<stk::mesh::EntityId>(chosenSideId);
}

typedef std::map<stk::mesh::GraphEdge, stk::mesh::impl::ParallelInfo *, GraphEdgeLessByElem2> GraphEdgeToParInfoMap;

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
    stk::mesh::impl::ParallelInfo& parInfoToChange = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(recvGraphEdge);
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

int get_elem2_proc_id(const stk::mesh::GraphEdge& graphEdge,
                      const stk::mesh::impl::ElementSidePair& MCEAndSide,
                      stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    int elem2Proc = -1;
    if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
    {
        stk::mesh::GraphEdge MCEtoElem2(MCEAndSide.first, MCEAndSide.second, graphEdge.elem2, graphEdge.side2);
        stk::mesh::impl::ParallelInfo& MCEtoElem2ParInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(MCEtoElem2);
        elem2Proc = MCEtoElem2ParInfo.get_proc_rank_of_neighbor();
    }
    return elem2Proc;
}

void pack_changed_edges_and_chosen_side_ids(stk::CommSparse& commSparse,
                                            const ElemSideAndEdges& elemSidesAndEdges,
                                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                            const IdMapper& idMapper)
{
    for(const ElemSideAndEdges::value_type& MCEsideAndEdges : elemSidesAndEdges)
    {
        const stk::mesh::impl::ElementSidePair& MCEAndSide = MCEsideAndEdges.first;
        const std::vector<stk::mesh::GraphEdge>& edgesForThisSide = MCEsideAndEdges.second;
        for(const stk::mesh::GraphEdge& graphEdge : edgesForThisSide)
        {
            if(!stk::mesh::impl::is_local_element(graphEdge.elem1))
            {
                stk::mesh::GraphEdge MCEtoElem1(MCEAndSide.first, MCEAndSide.second, graphEdge.elem1, graphEdge.side1);
                stk::mesh::impl::ParallelInfo& MCEtoElem1ParInfo = parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(MCEtoElem1);
                int elem1Proc = MCEtoElem1ParInfo.get_proc_rank_of_neighbor();
                int elem2Proc = get_elem2_proc_id(graphEdge, MCEAndSide, parallelInfoForGraphEdges);
                if(elem1Proc != elem2Proc)
                {
                    stk::mesh::GraphEdge edgeOnOtherProc(idMapper.localToGlobal(graphEdge.elem1),
                                                         graphEdge.side1,
                                                         idMapper.localToGlobal(graphEdge.elem2),
                                                         graphEdge.side2);
                    pack_graph_edge_and_chosen_side_id_to_proc(commSparse, elem1Proc, edgeOnOtherProc, MCEtoElem1ParInfo.m_chosen_side_id);
                }
            }
        }
    }
}

void update_chosen_ids_on_other_procs_for_edges_with_coincident_elements(const ElemSideAndEdges &elemSidesAndEdges,
                                        const stk::mesh::Graph &graph,
                                        const stk::mesh::impl::SparseGraph &coincidentEdges,
                                        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                        const IdMapper& idMapper,
                                        MPI_Comm comm)
{
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&commSparse, &elemSidesAndEdges, &parallelInfoForGraphEdges, &idMapper]()
        {
            pack_changed_edges_and_chosen_side_ids(commSparse, elemSidesAndEdges, parallelInfoForGraphEdges, idMapper);
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
    stk::mesh::impl::ParallelInfo& parInfoThisElem =
            parallelInfoForGraphEdges.get_parallel_info_for_graph_edge(parallelGraphEdge);
    parInfos.insert(std::make_pair(parallelGraphEdge, &parInfoThisElem));
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



void fill_par_infos_for_edges(const std::vector<stk::mesh::GraphEdge>& edges, stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges, GraphEdgeToParInfoMap &coincidentParInfos)
{
    for(const stk::mesh::GraphEdge& edge : edges)
        if(stk::mesh::impl::is_local_element(edge.elem1) && !stk::mesh::impl::is_local_element(edge.elem2))
            add_graph_edge_and_par_info(edge, parallelInfoForGraphEdges, coincidentParInfos);
}

std::vector<stk::mesh::GraphEdge> getMCEtoANCEedges(const stk::mesh::Graph& graph,
                                                    const stk::mesh::impl::LocalId MCE,
                                                    int sideIndex)
{
    std::vector<stk::mesh::GraphEdge> MCEtoANCEedges;
    for(const stk::mesh::GraphEdge& MCEtoANCEedge : graph.get_edges_for_element(MCE))
        if(MCEtoANCEedge.side1 == sideIndex)
            MCEtoANCEedges.push_back(MCEtoANCEedge);
    return MCEtoANCEedges;
}

std::vector<stk::mesh::GraphEdge> getMCEtoSCEedges(const std::vector<stk::mesh::GraphEdge>& allMCEtoSCEedges, int sideIndex)
{
    std::vector<stk::mesh::GraphEdge> MCEtoSCEedges;
    for(const stk::mesh::GraphEdge& MCEtoSCEedge : allMCEtoSCEedges)
        if(MCEtoSCEedge.side1 == sideIndex)
            MCEtoSCEedges.push_back(MCEtoSCEedge);
    return MCEtoSCEedges;
}

std::vector<stk::mesh::GraphEdge> getANCEtoSCEedges(const std::vector<stk::mesh::GraphEdge>& MCEtoANCEedges,
                                                    const std::vector<stk::mesh::GraphEdge>& allMCEtoSCEedges,
                                                    int sideIndex)
{
    std::vector<stk::mesh::GraphEdge> MCEtoSCEedges = getMCEtoSCEedges(allMCEtoSCEedges, sideIndex);

    std::vector<stk::mesh::GraphEdge> ANCEtoSCEedges;
    for(const stk::mesh::GraphEdge& MCEtoANCE : MCEtoANCEedges)
        for(const stk::mesh::GraphEdge& MCEtoSCE : MCEtoSCEedges)
            ANCEtoSCEedges.push_back(stk::mesh::GraphEdge(MCEtoANCE.elem2, MCEtoANCE.side2, MCEtoSCE.elem2, MCEtoSCE.side2));
    return ANCEtoSCEedges;
}

void appendTransposeEdges(const std::vector<stk::mesh::GraphEdge>& graphEdges, std::vector<stk::mesh::GraphEdge> &collectedEdges)
{
    for(const stk::mesh::GraphEdge& graphEdge : graphEdges)
        collectedEdges.push_back(stk::mesh::GraphEdge(graphEdge.elem2, graphEdge.side2, graphEdge.elem1, graphEdge.side1));
}

int get_num_sides_of_coincident_element(const std::vector<stk::mesh::GraphEdge>& allMCEtoSCEedges)
{
    int maxSideId = -1;
    for(const stk::mesh::GraphEdge& MCEtoSCEedge : allMCEtoSCEedges)
        maxSideId = std::max(maxSideId, MCEtoSCEedge.side1);
    return maxSideId + 1;
}

std::vector<stk::mesh::GraphEdge> get_all_edges_for_this_side(const stk::mesh::Graph& graph, const stk::mesh::impl::LocalId MCE , int sideIndex, const std::vector<stk::mesh::GraphEdge>& allMCEtoSCEedges)
{
    std::vector<stk::mesh::GraphEdge> allEdgesForThisSide;
    std::vector<stk::mesh::GraphEdge> MCEtoANCEedges = getMCEtoANCEedges(graph, MCE, sideIndex);
    std::vector<stk::mesh::GraphEdge> ANCEtoSCEedges = getANCEtoSCEedges(MCEtoANCEedges, allMCEtoSCEedges, sideIndex);

    allEdgesForThisSide.insert(allEdgesForThisSide.end(), MCEtoANCEedges.begin(), MCEtoANCEedges.end());
    allEdgesForThisSide.insert(allEdgesForThisSide.end(), ANCEtoSCEedges.begin(), ANCEtoSCEedges.end());
    appendTransposeEdges(MCEtoANCEedges, allEdgesForThisSide);
    appendTransposeEdges(ANCEtoSCEedges, allEdgesForThisSide);

    return allEdgesForThisSide;
}


void match_chosen_ids_for_edges_this_proc(const stk::mesh::Graph& graph,
                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                            const std::vector<stk::mesh::GraphEdge>& allMCEtoSCEedges,
                            const stk::mesh::impl::LocalId MCE,
                            const IdMapper &idMapper,
                            ElemSideAndEdges& elemSidesAndEdges,
                            MPI_Comm comm)
{
    int numMCEsides = get_num_sides_of_coincident_element(allMCEtoSCEedges);

    for(int sideIndex = 0; sideIndex < numMCEsides; sideIndex++)
    {
        std::vector<stk::mesh::GraphEdge> allEdgesForThisSide = get_all_edges_for_this_side(graph, MCE, sideIndex, allMCEtoSCEedges);

        GraphEdgeToParInfoMap coincidentParInfos;
        fill_par_infos_for_edges(allEdgesForThisSide, parallelInfoForGraphEdges, coincidentParInfos);

        stk::mesh::EntityId chosenSideId = get_min_chosen_side_id(coincidentParInfos);
        update_chosen_side_id_for_coincident_graph_edges(coincidentParInfos, chosenSideId);

        elemSidesAndEdges[stk::mesh::impl::ElementSidePair(MCE, sideIndex)] = allEdgesForThisSide;
    }
}

bool is_master_coincident_element(const stk::mesh::impl::LocalId elemId, const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem, const IdMapper &idMapper)
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
                                                 ElemSideAndEdges& elemSidesAndEdges,
                                                 MPI_Comm comm)
{
    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : extractedCoincidentElements)
    {
        const stk::mesh::impl::LocalId possibleMCE = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        if(is_master_coincident_element(possibleMCE, coincidentEdgesForElem, idMapper))
            match_chosen_ids_for_edges_this_proc(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, possibleMCE, idMapper, elemSidesAndEdges, comm);
    }
}

void make_chosen_ids_in_parinfo_consistent_for_edges_with_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parInfosForEdges,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm)
{
    ElemSideAndEdges elemSidesAndEdges;
    choose_consistent_face_ids_on_procs_that_own_coincident_elements(graph, parInfosForEdges, coincidentEdges, idMapper, elemSidesAndEdges, comm);
    update_chosen_ids_on_other_procs_for_edges_with_coincident_elements(elemSidesAndEdges, graph, coincidentEdges, parInfosForEdges, idMapper, comm);
}

}}} // end namespaces stk mesh impl


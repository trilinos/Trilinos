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


typedef std::map<stk::mesh::GraphEdge, stk::mesh::impl::parallel_info *, GraphEdgeLessByElem2> GraphEdgeToParInfoMap;

void pack_graph_edge_and_chosen_side_id(stk::CommSparse& commSparse,
                                        const stk::mesh::GraphEdge& graphEdge,
                                        const parallel_info* parInfo,
                                        const IdMapper& idMapper)
{
    stk::mesh::EntityId elem1GlobalId = idMapper.localToGlobal(graphEdge.elem1);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(-graphEdge.elem2);
    commSparse.send_buffer(parInfo->m_other_proc).pack<int>(graphEdge.side2);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(elem1GlobalId);
    commSparse.send_buffer(parInfo->m_other_proc).pack<int>(graphEdge.side1);
    commSparse.send_buffer(parInfo->m_other_proc).pack<stk::mesh::EntityId>(parInfo->m_chosen_side_id);
}

void pack_changed_edges_and_chosen_side_ids(stk::CommSparse& commSparse,
                                            const GraphEdgeToParInfoMap& changedEdgesAndParInfos,
                                            const IdMapper& idMapper)
{
    for(const GraphEdgeToParInfoMap::value_type& changedEdgeAndParInfo : changedEdgesAndParInfos)
    {
        const stk::mesh::GraphEdge& graphEdge = changedEdgeAndParInfo.first;
        const parallel_info* parInfo = changedEdgeAndParInfo.second;
        pack_graph_edge_and_chosen_side_id(commSparse, graphEdge, parInfo, idMapper);
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

void fix_chosen_side_ids_on_other_procs(const GraphEdgeToParInfoMap &changedEdgesAndParInfos,
                                        stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                        const IdMapper& idMapper,
                                        MPI_Comm comm)
{
    stk::CommSparse commSparse(comm);
    stk::pack_and_communicate(commSparse,
        [&commSparse, &changedEdgesAndParInfos, &idMapper]()
        {
            pack_changed_edges_and_chosen_side_ids(commSparse, changedEdgesAndParInfos, idMapper);
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
    return parallelGraphEdge.elem1 == coincidentGraphEdge.elem1 && parallelGraphEdge.side1 == coincidentGraphEdge.side1;
}

void add_coincident_graph_edges_and_par_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                              const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                                              stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                              GraphEdgeToParInfoMap& parInfos)
{
    for(const stk::mesh::GraphEdge& coincidentGraphEdge : coincidentEdgesForElem)
    {
        if(are_graph_edges_for_same_side(parallelGraphEdge, coincidentGraphEdge) &&
                impl::is_local_element(coincidentGraphEdge.elem2))
        {
            stk::mesh::GraphEdge otherElemGraphEdge(coincidentGraphEdge.elem2,
                                                    coincidentGraphEdge.side2,
                                                    parallelGraphEdge.elem2,
                                                    parallelGraphEdge.side2);
            add_graph_edge_and_par_info(otherElemGraphEdge, parallelInfoForGraphEdges, parInfos);
        }
    }
}

GraphEdgeToParInfoMap get_coincident_parallel_infos(const stk::mesh::GraphEdge& parallelGraphEdge,
                                                    const std::vector<stk::mesh::GraphEdge> &coincidentEdgesForElem,
                                                    stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges)
{
    GraphEdgeToParInfoMap parInfos;
    add_graph_edge_and_par_info(parallelGraphEdge, parallelInfoForGraphEdges, parInfos);
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
                                                      stk::mesh::EntityId chosenSideId,
                                                      GraphEdgeToParInfoMap& changedEdgesAndProcs)
{
    for(const GraphEdgeToParInfoMap::value_type &graphEdgeAndParInfo : coincidentParInfos)
    {
        if(graphEdgeAndParInfo.second->m_chosen_side_id != chosenSideId)
        {
            changedEdgesAndProcs.insert(graphEdgeAndParInfo);
            graphEdgeAndParInfo.second->m_chosen_side_id = chosenSideId;
        }
    }
}

void update_chosen_side_ids_for_remote_graph_edge(const stk::mesh::GraphEdge& graphEdge,
                                                  stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                  const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                                                  GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    GraphEdgeToParInfoMap coincidentParInfos = get_coincident_parallel_infos(graphEdge,
                                                                             coincidentEdgesForElem,
                                                                             parallelInfoForGraphEdges);
    stk::mesh::EntityId chosenSideId = get_min_chosen_side_id(coincidentParInfos);
    update_chosen_side_id_for_coincident_graph_edges(coincidentParInfos, chosenSideId, changedEdgesAndParInfos);
}

void update_chosen_side_ids(const stk::mesh::Graph& graph,
                            stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                            const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem,
                            const stk::mesh::impl::LocalId elemId,
                            GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    for(const stk::mesh::GraphEdge& graphEdge : graph.get_edges_for_element(elemId))
    {
        if(!stk::mesh::impl::is_local_element(graphEdge.elem2))
        {
            update_chosen_side_ids_for_remote_graph_edge(graphEdge,
                                                         parallelInfoForGraphEdges,
                                                         coincidentEdgesForElem,
                                                         changedEdgesAndParInfos);
        }
    }
}

void update_chosen_side_ids_for_coincident_elems(const stk::mesh::Graph& graph,
                                                 stk::mesh::ParallelInfoForGraphEdges& parallelInfoForGraphEdges,
                                                 const stk::mesh::impl::SparseGraph& extractedCoincidentElements,
                                                 GraphEdgeToParInfoMap& changedEdgesAndParInfos)
{
    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : extractedCoincidentElements)
    {
        const stk::mesh::impl::LocalId elemId = extractedEdgesForElem.first;
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        update_chosen_side_ids(graph, parallelInfoForGraphEdges, coincidentEdgesForElem, elemId, changedEdgesAndParInfos);
    }
}

void choose_face_id_for_coincident_elements(const stk::mesh::Graph &graph,
                                            stk::mesh::ParallelInfoForGraphEdges &parInfosForEdges,
                                            const stk::mesh::impl::SparseGraph &coincidentEdges,
                                            const IdMapper &idMapper,
                                            MPI_Comm comm)
{
    GraphEdgeToParInfoMap changedEdgesAndParInfos;
    update_chosen_side_ids_for_coincident_elems(graph, parInfosForEdges, coincidentEdges, changedEdgesAndParInfos);
    fix_chosen_side_ids_on_other_procs(changedEdgesAndParInfos, parInfosForEdges, idMapper, comm);
}

}}} // end namespaces stk mesh impl


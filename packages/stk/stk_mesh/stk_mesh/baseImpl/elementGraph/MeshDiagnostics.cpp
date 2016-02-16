#include <stddef.h>                     // for size_t, nullptr
#include <string>                       // for string
#include <stk_util/parallel/ParallelReduce.hpp>
#include <map>
#include <string>
#include "../../base/BulkData.hpp"
#include "../../base/MetaData.hpp"
#include "ElemElemGraph.hpp"
#include "MeshDiagnosticObserver.hpp"
#include "../EquivalentEntityBlocks.hpp"
#include "../../../../stk_util/stk_util/parallel/DistributedIndex.hpp"
#include "../../base/GetEntities.hpp"
#include "../MeshImplUtils.hpp"

namespace stk { namespace mesh {

class ElemGraphForDiagnostics : public stk::mesh::ElemElemGraph
{
public:
    ElemGraphForDiagnostics(stk::mesh::BulkData& bulkData, const stk::mesh::Selector &selector, const stk::mesh::Selector *air = nullptr) :
        stk::mesh::ElemElemGraph(bulkData, selector, air) { }

    stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
    stk::mesh::Entity get_entity(stk::mesh::impl::LocalId local_id) const { return m_local_id_to_element_entity[local_id]; }
    const stk::mesh::impl::SparseGraph& get_coincident_graph() const { return m_coincidentGraph; }
};

std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > get_split_coincident_elements(stk::mesh::BulkData& bulkData)
{
    stk::mesh::Selector sel = bulkData.mesh_meta_data().locally_owned_part();
    ElemGraphForDiagnostics graph(bulkData, sel);
    const stk::mesh::impl::SparseGraph& coingraph = graph.get_coincident_graph();

    std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > badElements;

    for(const stk::mesh::impl::SparseGraph::value_type& extractedEdgesForElem : coingraph)
    {
        const std::vector<stk::mesh::GraphEdge>& coincidentEdgesForElem = extractedEdgesForElem.second;
        for(const stk::mesh::GraphEdge& edge : coincidentEdgesForElem)
        {
            if(edge.elem2 < 0)
            {
                stk::mesh::Entity entity = graph.get_entity(edge.elem1);
                stk::mesh::EntityId id = bulkData.identifier(entity);
                stk::mesh::impl::ParallelGraphInfo& par_info = graph.get_parallel_info();
                stk::mesh::impl::ParallelGraphInfo::iterator iter = par_info.find(edge);
                ThrowRequireMsg(iter!=par_info.end(), "Program error. Contact sierra-help@sandia.gov for support.");
                badElements[id] = std::make_pair(-edge.elem2, iter->second.m_other_proc);
            }
        }
    }
    return badElements;
}

std::string get_message_for_split_coincident_elements(const stk::mesh::BulkData& bulkData, const std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > & splitCoincidentElements)
{
    std::ostringstream out;
    for(const auto& item : splitCoincidentElements) {
        stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK,item.first);
        const stk::mesh::PartVector& elementParts = bulkData.bucket(element).supersets();
        std::string blockNames;
        blockNames = "{";
        for (const stk::mesh::Part* part : elementParts) {
            if (stk::mesh::impl::is_element_block(*part)) {
                blockNames += " " + part->name();
            }
        }
        blockNames += " }";
        out << "[" << bulkData.parallel_rank() << "] Element " << item.first << " (" << bulkData.bucket(element).topology() << ") in blocks " << blockNames << " is coincident with element " << item.second.first << " on processor " << item.second.second << std::endl;
    }
    return out.str();
}

void throw_if_any_proc_has_false(MPI_Comm comm, bool is_all_ok_locally)
{
    bool is_all_ok_globally = stk::is_true_on_all_procs(comm, is_all_ok_locally);
    ThrowRequireMsg(is_all_ok_globally, "Mesh diagnostics failed.");
}

stk::mesh::Selector get_owned_or_shared_selector(const stk::mesh::BulkData & bulkData)
{
    return bulkData.mesh_meta_data().locally_owned_part() | bulkData.mesh_meta_data().globally_shared_part();
}

stk::parallel::DistributedIndex::KeyTypeVector get_all_local_keys(const stk::mesh::BulkData & bulkData)
{
    stk::parallel::DistributedIndex::KeyTypeVector localKeys;
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK;rank < bulkData.mesh_meta_data().entity_rank_count();++rank)
    {
        stk::mesh::EntityVector entities;
        stk::mesh::get_selected_entities(get_owned_or_shared_selector(bulkData), bulkData.buckets(rank), entities);
        for(stk::mesh::Entity entity: entities)
            localKeys.push_back(bulkData.entity_key(entity));
    }
    return localKeys;
}

void add_keys_to_distributed_index(const stk::mesh::BulkData & bulkData, stk::parallel::DistributedIndex & distributedIndex)
{
    stk::parallel::DistributedIndex::KeyTypeVector localKeys = get_all_local_keys(bulkData);

    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator begin = localKeys.begin();
    stk::parallel::DistributedIndex::KeyTypeVector::const_iterator end = localKeys.end();
    distributedIndex.update_keys( begin, end );
}

std::vector<stk::mesh::EntityKeyProc> get_non_unique_keys(const stk::mesh::BulkData& bulkData, const stk::parallel::DistributedIndex& distributedIndex,
        const stk::parallel::DistributedIndex::KeyTypeVector& localKeys)
{
    stk::parallel::DistributedIndex::KeyProcVector sharedKeyProcs;
    distributedIndex.query_to_usage(localKeys, sharedKeyProcs);

    std::vector<stk::mesh::EntityKeyProc> badKeys;
    for (const stk::parallel::DistributedIndex::KeyProc& sharedKeyProc : sharedKeyProcs)
    {
        stk::mesh::EntityKey key( static_cast<stk::mesh::EntityKey::entity_key_t>(sharedKeyProc.first) );
        if ( bulkData.parallel_rank() != sharedKeyProc.second )
        {
            if(!bulkData.in_shared(key, sharedKeyProc.second))
                badKeys.push_back({key, sharedKeyProc.second});
        }
    }
    return badKeys;
}

std::string get_topology(stk::topology topology)
{
    if(topology==stk::topology::INVALID_TOPOLOGY)
        return " ";
    return " (" + topology.name() + ") ";
}

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData)
{
    stk::parallel::DistributedIndex distributedIndex( bulkData.parallel(), stk::mesh::impl::convert_entity_keys_to_spans(bulkData.mesh_meta_data()));
    add_keys_to_distributed_index(bulkData, distributedIndex);
    stk::parallel::DistributedIndex::KeyTypeVector localKeys = get_all_local_keys(bulkData);
    return get_non_unique_keys(bulkData, distributedIndex, localKeys);
}

std::string get_non_unique_key_messages(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs)
{
    std::ostringstream os;
    for(const stk::mesh::EntityKeyProc& keyProc : badKeyProcs)
    {
        stk::mesh::Entity entity = bulkData.get_entity(keyProc.first);
        os << "[" << bulkData.parallel_rank() << "] Key " << keyProc.first <<
                get_topology(bulkData.bucket(entity).topology()) << "is also present (inappropriately) on processor " <<
                keyProc.second << "." << std::endl;
    }
    return os.str();
}


} }

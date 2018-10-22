#include <cmath>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/GeometricVertices.hpp>
#include <stk_balance/internal/StkGeometricMethodViaZoltan.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/FieldBase.hpp>  // for field_data
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <zoltan.h>
#include <Zoltan2_Version.hpp>
#include <map>

#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include <stk_util/environment/Env.hpp>
#include "stk_mesh/base/FieldParallel.hpp"

namespace stk {
namespace balance {

class DecompositionChangeList::Impl
{
public:
    Impl(stk::mesh::BulkData & stkMeshBulkData, const stk::mesh::EntityProcVec & decomposition)
      : m_stkMeshBulkData(stkMeshBulkData)
    {
        fill_change_list_from_raw_decomposition(decomposition);
    }

    bool has_entity(stk::mesh::Entity entity) const
    {
        return (m_dataMap.find(entity) != m_dataMap.end());
    }

    int get_entity_destination(stk::mesh::Entity entity) const
    {
        std::map<stk::mesh::Entity, int>::const_iterator entry = m_dataMap.find(entity);
        if (entry != m_dataMap.end()) {
            return entry->second;
        }
        return -1;
    }

    void set_entity_destination(stk::mesh::Entity entity, const int destination)
    {
        m_dataMap[entity] = destination;
    }

    void delete_entity(stk::mesh::Entity entity)
    {
        std::map<stk::mesh::Entity, int>::iterator entry = m_dataMap.find(entity);

        if (entry != m_dataMap.end()) {
            m_dataMap.erase(entry);
        }
    }

    stk::mesh::EntityProcVec get_all_partition_changes()
    {
        stk::mesh::EntityProcVec entityProcPairs;
        for (auto & entry : m_dataMap) {
            entityProcPairs.push_back(entry);
        }
        return entityProcPairs;
    }

    stk::mesh::BulkData & get_bulk() { return m_stkMeshBulkData; }

    void get_decomposition_with_full_closure(stk::mesh::EntityProcVec & finalDecomposition)
    {
        finalDecomposition.clear();
        for (auto & entry : m_dataMap) {
            finalDecomposition.push_back(entry);
            add_downward_closure_for_entity(entry, finalDecomposition);
        }
    }

    size_t get_num_global_entity_migrations() const
    {
        size_t num_local_entity_migrations = m_dataMap.size();
        size_t num_global_entity_migrations = 0;
        stk::all_reduce_sum(m_stkMeshBulkData.parallel(), &num_local_entity_migrations, &num_global_entity_migrations, 1);
        return num_global_entity_migrations;
    }

    size_t get_max_global_entity_migrations() const
    {
        size_t num_local_entity_migrations = m_dataMap.size();
        size_t max_global_entity_migrations = 0;
        stk::all_reduce_max(m_stkMeshBulkData.parallel(), &num_local_entity_migrations, &max_global_entity_migrations, 1);
        return max_global_entity_migrations;
    }

private:
    stk::mesh::BulkData &m_stkMeshBulkData;
    std::map<stk::mesh::Entity, int> m_dataMap;

    void fill_change_list_from_raw_decomposition(const stk::mesh::EntityProcVec& decomposition) {
        for(const stk::mesh::EntityProc& entity_proc : decomposition)
        {
            if(m_stkMeshBulkData.is_valid(entity_proc.first) && entity_proc.second != m_stkMeshBulkData.parallel_rank())
                m_dataMap[entity_proc.first] = entity_proc.second;
        }
    }

    void add_downward_closure_for_entity(const std::pair<const stk::mesh::Entity, int> & entityProc, stk::mesh::EntityProcVec & finalDecomposition)
    {
        const stk::topology::rank_t entityRank = m_stkMeshBulkData.entity_rank(entityProc.first);
        for (int rank = entityRank-1; rank >= stk::topology::NODE_RANK; --rank)
            internal::add_connected_entities_of_rank(m_stkMeshBulkData, entityProc.first, entityProc.second, static_cast<stk::topology::rank_t>(rank), finalDecomposition);
    }

};


DecompositionChangeList::DecompositionChangeList(stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::EntityProcVec& decomposition)
    : pImpl(new Impl(stkMeshBulkData, decomposition))
{ }
DecompositionChangeList::~DecompositionChangeList() { delete pImpl;  }

bool DecompositionChangeList::has_entity(stk::mesh::Entity entity)                                    { return pImpl->has_entity(entity); }
int  DecompositionChangeList::get_entity_destination(stk::mesh::Entity entity)                        { return pImpl->get_entity_destination(entity); }
void DecompositionChangeList::set_entity_destination(stk::mesh::Entity entity, const int destination) { pImpl->set_entity_destination(entity, destination); }
void DecompositionChangeList::delete_entity(stk::mesh::Entity entity)                                 { pImpl->delete_entity(entity); }
stk::mesh::BulkData & DecompositionChangeList::get_bulk()                                             { return pImpl->get_bulk(); }
stk::mesh::EntityProcVec DecompositionChangeList::get_all_partition_changes()                         { return pImpl->get_all_partition_changes(); }
size_t DecompositionChangeList::get_num_global_entity_migrations() const                              { return pImpl->get_num_global_entity_migrations(); }
size_t DecompositionChangeList::get_max_global_entity_migrations() const                              { return pImpl->get_max_global_entity_migrations(); }

namespace internal
{


unsigned get_local_id(const stk::mesh::impl::LocalIdMapper& localIds, stk::mesh::Entity entity)
{
    return localIds.entity_to_local(entity);
}

void addBoxForFace(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity face, const double eps, BoxVectorWithStkId &faceBoxes, const stk::mesh::FieldBase* coord)
{

    unsigned numElements = stkMeshBulkData.num_elements(face);
    const stk::mesh::Entity *element = stkMeshBulkData.begin_elements(face);

    ThrowRequireWithSierraHelpMsg(numElements <= 1);
    if ( element != NULL && stkMeshBulkData.bucket(*element).owned() )
    {
        unsigned numNodes = stkMeshBulkData.num_nodes(face);
        const stk::mesh::Entity *nodes = stkMeshBulkData.begin_nodes(face);

        addBoxForNodes(stkMeshBulkData, numNodes, nodes, coord, eps, stkMeshBulkData.identifier(*element), faceBoxes);
    }
}

void addEdgeAndVertexWeightsForSearchResult(stk::mesh::BulkData& stkMeshBulkData, const BalanceSettings &balanceSettings, stk::mesh::EntityId element1Id,
        stk::mesh::EntityId element2Id, unsigned owningProcElement2, std::vector<GraphEdge>& graphEdges)
{
    stk::mesh::EntityKey entityKeyElement1(stk::topology::ELEMENT_RANK, element1Id);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(entityKeyElement1);
    ThrowRequireWithSierraHelpMsg(stkMeshBulkData.entity_rank(element1) == stk::topology::ELEMENT_RANK);

    if(stkMeshBulkData.is_valid(element1) && stkMeshBulkData.bucket(element1).owned() && element1Id != element2Id)
    {
        stk::mesh::EntityKey entityKeyElement2(stk::topology::ELEMENT_RANK, element2Id);
        stk::mesh::Entity element2 = stkMeshBulkData.get_entity(entityKeyElement2);

        unsigned anyIntersections = 0;
        if ( stkMeshBulkData.is_valid(element2) )
        {
            anyIntersections = stk::balance::internal::getNumSharedNodesBetweenElements(stkMeshBulkData, element1, element2);
        }

        bool elementIsNotConnectedViaNodes = anyIntersections == 0;
        if ( elementIsNotConnectedViaNodes )
        {
            double edge_weight = balanceSettings.getGraphEdgeWeightForSearch();
            graphEdges.push_back(GraphEdge(element1, element2Id, owningProcElement2, edge_weight, true));
        }
    }
}

void createGraphEdgesUsingBBSearch(stk::mesh::BulkData& stkMeshBulkData, const BalanceSettings &balanceSettings, std::vector<GraphEdge>& graphEdges,
                                   const stk::mesh::Selector& searchSelector)
{
    std::ostringstream os;
    size_t max = 0, min = 0, avg = 0;
    stk::get_max_min_avg(stkMeshBulkData.parallel(), graphEdges.size(), max, min, avg);
    os << "Starting search, have following distribution of graph edges: min=" << min << ", avg=" << avg << ", max=" << max;
    logMessage(stkMeshBulkData.parallel(), os.str());
    os.str("");

    StkSearchResults searchResults = stk::balance::internal::getSearchResultsForFacesParticles(stkMeshBulkData, balanceSettings, searchSelector);

    stk::get_max_min_avg(stkMeshBulkData.parallel(), searchResults.size(), max, min, avg);
    os << "Finished search, have following distribution of search results: min=" << min << ", avg=" << avg << ", max=" << max;
    logMessage(stkMeshBulkData.parallel(), os.str());
    os.str("");

    for(size_t i = 0; i < searchResults.size(); i++)
    {
        stk::mesh::EntityId element1Id = searchResults[i].first.id();
        stk::mesh::EntityId element2Id = searchResults[i].second.id();
        int owningProcElement1 = searchResults[i].first.proc();
        int owningProcElement2 = searchResults[i].second.proc();

        ThrowRequireWithSierraHelpMsg(owningProcElement1 == stkMeshBulkData.parallel_rank() || owningProcElement2 == stkMeshBulkData.parallel_rank());

        if ( owningProcElement1 == stkMeshBulkData.parallel_rank() )
        {
            addEdgeAndVertexWeightsForSearchResult(stkMeshBulkData, balanceSettings, element1Id, element2Id, owningProcElement2, graphEdges);
        }
        else
        {
            addEdgeAndVertexWeightsForSearchResult(stkMeshBulkData, balanceSettings, element2Id, element1Id, owningProcElement1, graphEdges);
        }
    }

    stk::get_max_min_avg(stkMeshBulkData.parallel(), graphEdges.size(), max, min, avg);
    os << "After search, have following distribution of graph edges: min=" << min << ", avg=" << avg << ", max=" << max;
    logMessage(stkMeshBulkData.parallel(), os.str());
}

std::vector<int> getLocalIdsOfEntitiesNotSelected(const stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Selector selector, const stk::mesh::impl::LocalIdMapper& localIds)
{
    selector = (!selector) & stkMeshBulkData.mesh_meta_data().locally_owned_part();
    std::vector<int> local_ids;
    size_t num_total_elements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

    const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, selector);
    for (size_t i=0; i<buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        for (size_t j=0;j<bucket.size();j++)
        {
            unsigned local_id = get_local_id(localIds, bucket[j]);
            ThrowRequireWithSierraHelpMsg(local_id < num_total_elements);
            local_ids.push_back(local_id);
        }
    }
    return local_ids;
}

void fillEntityCentroid(const stk::mesh::BulkData &stkMeshBulkData, const stk::mesh::FieldBase* coord, stk::mesh::Entity entityOfConcern, double *entityCentroid)
{
    unsigned spatialDimension = stkMeshBulkData.mesh_meta_data().spatial_dimension();
    if(stkMeshBulkData.entity_rank(entityOfConcern)==stk::topology::NODE_RANK)
    {
        double *nodeCoord = static_cast<double *>(stk::mesh::field_data(*coord, entityOfConcern));
        for(unsigned k=0; k < spatialDimension; k++)
        {
            entityCentroid[k] = nodeCoord[k];
        }
    }
    else
    {
        stk::mesh::Entity const * nodes = stkMeshBulkData.begin_nodes(entityOfConcern);
        unsigned numNodes = stkMeshBulkData.num_nodes(entityOfConcern);
        for(unsigned nodeIndex=0; nodeIndex < numNodes; nodeIndex++)
        {
            double *nodeCoord = static_cast<double *>(stk::mesh::field_data(*coord, nodes[nodeIndex]));
            for(unsigned k=0; k < spatialDimension; k++)
            {
                entityCentroid[k] += nodeCoord[k];
            }
        }
        for(unsigned k=0; k < spatialDimension; k++)
        {
            entityCentroid[k] /= numNodes;
        }
    }
}

void fill_list_of_entities_to_send_for_aura_like_ghosting(stk::mesh::BulkData& bulkData, stk::mesh::EntityProcVec &entitiesToGhost)
{
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part());
    for(const stk::mesh::Bucket *bucket : buckets)
    {
        for(stk::mesh::Entity shared_node : *bucket)
        {
            unsigned num_elements = bulkData.num_elements(shared_node);
            const stk::mesh::Entity* elements = bulkData.begin_elements(shared_node);
            for(unsigned i=0;i<num_elements;++i)
            {
                if(bulkData.bucket(elements[i]).owned())
                {
                    std::vector<int> comm_shared_procs;
                    bulkData.comm_shared_procs(bulkData.entity_key(shared_node), comm_shared_procs);
                    for(int proc : comm_shared_procs )
                    {
                        entitiesToGhost.push_back(stk::mesh::EntityProc(elements[i], proc));
                    }
                }
            }
        }
    }
}

void fill_connectivity_count_field(stk::mesh::BulkData & stkMeshBulkData, const BalanceSettings & balanceSettings)
{
    if (balanceSettings.shouldFixSpiders()) {
        const stk::mesh::Field<int> * connectivityCountField = balanceSettings.getSpiderConnectivityCountField(stkMeshBulkData);

        stk::mesh::Selector beamNodesSelector(stkMeshBulkData.mesh_meta_data().locally_owned_part() &
                                              stkMeshBulkData.mesh_meta_data().get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::BEAM_2)));
        const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, beamNodesSelector);

        for (stk::mesh::Bucket * bucket : buckets) {
            for (stk::mesh::Entity node : *bucket) {
                int * connectivityCount = stk::mesh::field_data(*connectivityCountField, node);
                const unsigned numElements = stkMeshBulkData.num_elements(node);
                const stk::mesh::Entity *element = stkMeshBulkData.begin_elements(node);
                for (unsigned i = 0; i < numElements; ++i) {
                    if (stkMeshBulkData.bucket(element[i]).topology() == stk::topology::BEAM_2) {
                        (*connectivityCount)++;
                    }
                }
            }
        }

        stk::mesh::communicate_field_data(stkMeshBulkData, {connectivityCountField});
    }
}

stk::mesh::Ghosting * create_custom_ghosting(stk::mesh::BulkData & stkMeshBulkData, const BalanceSettings & balanceSettings)
{
    stk::mesh::Ghosting * customAura = nullptr;
    if (!stkMeshBulkData.is_automatic_aura_on())
    {
        stkMeshBulkData.modification_begin();

        customAura = &stkMeshBulkData.create_ghosting("customAura");
        stk::mesh::EntityProcVec entitiesToGhost;
        fill_list_of_entities_to_send_for_aura_like_ghosting(stkMeshBulkData, entitiesToGhost);
        stkMeshBulkData.change_ghosting(*customAura, entitiesToGhost);

        stkMeshBulkData.modification_end();
    }

    return customAura;
}

void destroy_custom_ghosting(stk::mesh::BulkData & stkMeshBulkData, stk::mesh::Ghosting * customAura)
{
    if (nullptr != customAura)
    {
        stkMeshBulkData.modification_begin();
        stkMeshBulkData.destroy_ghosting(*customAura);
        stkMeshBulkData.modification_end();
    }
}

void logMessage(MPI_Comm communicator, const std::string &message)
{
    static double startTime = stk::wall_time();
    double now = stk::wall_time();

    size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;
    stk::get_memory_high_water_mark_across_processors(communicator, hwm_max, hwm_min, hwm_avg);

    sierra::Env::outputP0() << "[time:" << std::fixed << std::setprecision(3) << std::setw(10) << now-startTime << " s, hwm:"
              << std::setfill(' ') << std::right << std::setw(8) << stk::human_bytes(hwm_avg) << "] "<< message << std::endl;
}



void fill_zoltan2_graph(const BalanceSettings& balanceSettings,
                        stk::mesh::BulkData& stkMeshBulkData,
                        Zoltan2ParallelGraph& zoltan2Graph,
                        const stk::mesh::Selector& searchSelector,
                        const stk::mesh::impl::LocalIdMapper& localIds)
{
    std::vector<int> adjacencyProcs;
    zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                   balanceSettings,
                                                   adjacencyProcs,
                                                   searchSelector,
                                                   localIds);

    logMessage(stkMeshBulkData.parallel(), "Finished filling in graph data");
}

bool is_geometric_method(const std::string method)
{
    return (method=="rcb" || method=="rib" || method=="multijagged");
}


//#define WRITE_OUT_DEBUGGING_INFO
//#define WRITE_OUT_DECOMP_METRICS#include "StkGeometricMethodViaZoltan.hpp"


stk::mesh::EntityVector get_entities_to_balance(stk::mesh::Selector selector, stk::mesh::EntityRank primaryRank, const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector entitiesToBalance;
    selector = selector & bulkData.mesh_meta_data().locally_owned_part();
    stk::mesh::get_selected_entities(selector, bulkData.buckets(primaryRank), entitiesToBalance);
    return entitiesToBalance;
}


void get_multicriteria_decomp_using_selectors_as_segregation(const stk::mesh::BulkData& mesh, const std::vector<stk::mesh::Selector>& criterions, const BalanceSettings& balanceSettings,
                                                             const int numSubdomainsToCreate, stk::mesh::EntityProcVec& decomp, const stk::mesh::impl::LocalIdMapper& localIds)
{
    stk::mesh::Selector unionSelector = stk::mesh::selectUnion(criterions);
    stk::mesh::EntityVector entitiesToBalance = get_entities_to_balance(unionSelector, stk::topology::ELEM_RANK, mesh);
    size_t num_entities = entitiesToBalance.size();
    size_t num_entities_across_procs = 0;
    stk::all_reduce_sum(mesh.parallel(), &num_entities, &num_entities_across_procs, 1);
    if(num_entities_across_procs > 0)
    {
        stk::balance::internal::GeometricVertices vertexInfo(balanceSettings, mesh, entitiesToBalance, criterions);
        std::vector<unsigned> processorOntoWhichEntityBelongs = stk::balance::get_decomposition(vertexInfo, balanceSettings, numSubdomainsToCreate, mesh.parallel());
        for(size_t i=0;i<entitiesToBalance.size();++i)
        {
            int local_id = get_local_id(localIds, entitiesToBalance[i]);
            decomp[local_id] = std::make_pair(entitiesToBalance[i], processorOntoWhichEntityBelongs[i]);
        }
    }
}


void fill_decomp_using_geometric_method(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate, stk::mesh::EntityProcVec &decomp,
                                        stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds)
{
    logMessage(stkMeshBulkData.parallel(), "Using Zoltan2 version: " + Zoltan2::Zoltan2_Version());
    logMessage(stkMeshBulkData.parallel(), "Filling in vertex data for decomp method = " + balanceSettings.getDecompMethod());

    size_t numEntities = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

    decomp.clear();
    decomp.resize(numEntities, std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));


    if (balanceSettings.isMultiCriteriaRebalance())
        get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, selectors, balanceSettings, numSubdomainsToCreate, decomp, localIds);
    else
        for(const stk::mesh::Selector selector : selectors)
            get_multicriteria_decomp_using_selectors_as_segregation(stkMeshBulkData, std::vector<stk::mesh::Selector>{selector}, balanceSettings, numSubdomainsToCreate, decomp, localIds);

    logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");
}

void get_multicriteria_parmetis_decomp(const stk::mesh::BulkData &mesh, const BalanceSettings& balanceSettings, Zoltan2ParallelGraph &zoltan2Graph, Teuchos::ParameterList &params,
                                       stk::mesh::Selector selector, stk::mesh::EntityProcVec &decomp, const stk::mesh::impl::LocalIdMapper& localIds)
{
    StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

    logMessage(mesh.parallel(), "Setting up partitioning problem");

    Zoltan2::PartitioningProblem<StkMeshZoltanAdapter> problem(&stkMeshAdapter, &params, mesh.parallel());

    logMessage(mesh.parallel(), "Solving");

    if(balanceSettings.shouldPrintMetrics())
    {
        internal::print_statistics(stkMeshAdapter, mesh.parallel(), mesh.parallel_rank());
    }

    std::srand(mesh.parallel_rank()); // KHP: Temporary until an API is added to Zoltan2 for random seeds.
    problem.solve();

    if(balanceSettings.shouldPrintMetrics())
        internal::print_solution_statistics(stkMeshAdapter, problem.getSolution(), mesh.parallel(), mesh.parallel_rank());

#if defined(WRITE_OUT_DEBUGGING_INFO)
    std::vector<int> local_ids_of_elements_to_balance;
    std::set_difference(all_local_ids.begin(), all_local_ids.end(), localIdsOfNotSelectedEntities.begin(), localIdsOfNotSelectedEntities.end(),
                            std::inserter(local_ids_of_elements_to_balance, local_ids_of_elements_to_balance.begin()));

    std::ostringstream filename;
    filename << "rebalance_proc_" << mesh.parallel_rank() << "_for_selector_" << i << "_for_step_" << step << ".txt";
    std::ofstream out(filename.str().c_str());
    out << "For this selector: " << selectors[i] << " the following entities are being balanced: ";
    for(size_t j=0;j<local_ids_of_elements_to_balance.size();++j)
    {
        stk::mesh::Entity element = entities[local_ids_of_elements_to_balance[j]];
        out << "Element " << j << " has: " << mesh.entity_key(element) << " with topology " << mesh.bucket(element).topology() << std::endl;
    }
    stkMeshAdapter.debuggingInfo(mesh.parallel_rank(), out);
    out.close();
#endif

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(mesh.mesh_meta_data().locally_owned_part() & selector, mesh.buckets(stk::topology::ELEM_RANK), elements);

    // local entity j is on which processor
    const StkMeshZoltanAdapter::part_t *processorOntoWhichEntityBelongs = problem.getSolution().getPartListView();
    for(size_t j = 0; j < elements.size(); ++j)
    {
        int local_id = get_local_id(localIds, elements[j]);
        int dest_proc = processorOntoWhichEntityBelongs[local_id];
        decomp[local_id] = std::make_pair(elements[j], dest_proc);
    }
}

Teuchos::ParameterList getGraphBasedParameters(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate)
{
    Teuchos::ParameterList params("test params");
    params.set("debug_level", "no_status");
//    params.set("debug_level", "basic_status");

    int nparts = numSubdomainsToCreate;
    double imbalance_allowed = balanceSettings.getImbalanceTolerance();
    params.set("imbalance_tolerance", imbalance_allowed);
    params.set("num_global_parts", nparts);
    params.set("algorithm", balanceSettings.getDecompMethod());

//    params.set("partitioning_objective", "balance_object_weight");
//    params.set("partitioning_objective", "multicriteria_minimize_total_weight");
//    params.set("partitioning_objective", "multicriteria_minimize_maximum_weight");
//    params.set("partitioning_objective", "multicriteria_balance_total_maximum");

    if (balanceSettings.isIncrementalRebalance())
    {
        params.set("partitioning_approach", "repartition");
        params.set("remap_parts", true);
    }

    // should not hurt other methods, only affects RCB.
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters", false);
    zparams.set("debug_level", "0");
//    zparams.set("LB_METHOD", "PHG");
//    zparams.set("LB_APPROACH", "PARTITION");
    return params;
}

bool isElementPartOfSpider(const stk::mesh::BulkData& stkMeshBulkData,
                           const stk::mesh::Field<int>& beamConnectivityCountField,
                           stk::mesh::Entity element)
{
    const int spiderConnectivityThreshold = 5;
    const stk::mesh::Entity* nodes = stkMeshBulkData.begin_nodes(element);
    const unsigned numNodes = stkMeshBulkData.num_nodes(element);
    for (unsigned i = 0; i < numNodes; ++i) {
        const int connectivityCount = *stk::mesh::field_data(beamConnectivityCountField, nodes[i]);
        if (connectivityCount > spiderConnectivityThreshold) {
            return true;
        }
    }
    return false;
}

bool shouldOmitSpiderElement(const stk::mesh::BulkData & stkMeshBulkData,
                             const stk::balance::BalanceSettings & balanceSettings,
                             stk::mesh::Entity elem)
{
    bool omitConnection = false;
    if (balanceSettings.shouldFixSpiders()) {
        const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderConnectivityCountField(stkMeshBulkData);
        stk::topology elemTopology = stkMeshBulkData.bucket(elem).topology();

        if (elemTopology == stk::topology::BEAM_2 || elemTopology == stk::topology::PARTICLE) {
            omitConnection = isElementPartOfSpider(stkMeshBulkData, beamConnectivityCountField, elem);
        }
    }

    return omitConnection;
}

void fix_spider_elements(const BalanceSettings & balanceSettings, stk::mesh::BulkData & stkMeshBulkData)
{
    stk::mesh::Ghosting * customAura = create_custom_ghosting(stkMeshBulkData, balanceSettings);

    stk::mesh::MetaData & meta = stkMeshBulkData.mesh_meta_data();
    const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderConnectivityCountField(stkMeshBulkData);

    stk::mesh::EntityVector beams;
    stk::mesh::Part & beamPart = meta.get_topology_root_part(stk::topology::BEAM_2);
    stk::mesh::get_selected_entities(beamPart & meta.locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK), beams);

    stk::mesh::EntityProcVec beamsToMove;
    for (stk::mesh::Entity beam : beams) {
        if (isElementPartOfSpider(stkMeshBulkData, beamConnectivityCountField, beam)) {
            const stk::mesh::Entity* nodes = stkMeshBulkData.begin_nodes(beam);
            const int node1ConnectivityCount = *stk::mesh::field_data(beamConnectivityCountField, nodes[0]);
            const int node2ConnectivityCount = *stk::mesh::field_data(beamConnectivityCountField, nodes[1]);

            const stk::mesh::Entity endNode = (node1ConnectivityCount < node2ConnectivityCount) ? nodes[0] : nodes[1];
            const stk::mesh::Entity* elements = stkMeshBulkData.begin_elements(endNode);
            const unsigned numElements = stkMeshBulkData.num_elements(endNode);
            int newOwner = stkMeshBulkData.parallel_size() - 1;
            for (unsigned i = 0; i < numElements; ++i) {
                if (stkMeshBulkData.bucket(elements[i]).topology() != stk::topology::BEAM_2) {
                    newOwner = std::min(newOwner, stkMeshBulkData.parallel_owner_rank(elements[i]));
                }
            }

            if (newOwner != stkMeshBulkData.parallel_rank()) {
                beamsToMove.push_back(std::make_pair(beam, newOwner));
            }
        }
    }

    destroy_custom_ghosting(stkMeshBulkData, customAura);

    stkMeshBulkData.change_entity_owner(beamsToMove);
}

void keep_spiders_on_original_proc(stk::mesh::BulkData &bulk, const stk::balance::BalanceSettings & balanceSettings, DecompositionChangeList &changeList)
{
    // Need to keep spiders on the original proc until the remaining elements have moved,
    // so that we can properly determine the final ownership of the elements on the end.
    // Then, we can move them.
    //
    const stk::mesh::Field<int> & beamConnectivityCountField = *balanceSettings.getSpiderConnectivityCountField(bulk);

    stk::mesh::EntityProcVec entityProcs = changeList.get_all_partition_changes();
    for (const stk::mesh::EntityProc & entityProc : entityProcs) {
        stk::mesh::Entity entity = entityProc.first;
        if (bulk.bucket(entity).topology() == stk::topology::BEAM_2) {
            if (isElementPartOfSpider(bulk, beamConnectivityCountField, entity)) {
                changeList.delete_entity(entity);
            }
        }
    }
}

void createZoltanParallelGraph(const BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData,
                               const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds,
                               Zoltan2ParallelGraph& zoltan2Graph)
{
    std::vector<size_t> counts;
    stk::mesh::Selector locallyOwnedSelector(stkMeshBulkData.mesh_meta_data().locally_owned_part());

    stk::mesh::comm_mesh_counts(stkMeshBulkData, counts, &locallyOwnedSelector);
    zoltan2Graph.set_num_global_elements(counts[stk::topology::ELEM_RANK]);
    zoltan2Graph.set_spatial_dim(stkMeshBulkData.mesh_meta_data().spatial_dimension());
    if (balanceSettings.isMultiCriteriaRebalance())
        zoltan2Graph.set_num_field_criteria(selectors.size()*balanceSettings.getNumCriteria());

    stk::mesh::Selector selectUnion = stk::mesh::selectUnion(selectors);

    // set vertex weights using entity's topology and if search is part of algorithm, use multiplier
    fill_zoltan2_graph(balanceSettings, stkMeshBulkData, zoltan2Graph, selectUnion, localIds);

    // now can reset those vertex weights based on fields or other critieria
    zoltan2Graph.adjust_vertex_weights(balanceSettings, stkMeshBulkData, selectors, localIds);

    if (balanceSettings.allowModificationOfVertexWeightsForSmallMeshes())
    {
        bool isSmallMesh = (counts[stk::topology::ELEM_RANK] / stkMeshBulkData.parallel_size()) <= 10;
        if(isSmallMesh)
        {
            logMessage(stkMeshBulkData.parallel(), "Changing weights since mesh is small");
            zoltan2Graph.adjust_weights_for_small_meshes();
        }
    }
}

void fill_decomp_using_parmetis(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate, stk::mesh::EntityProcVec &decomp, stk::mesh::BulkData& stkMeshBulkData,
                                const std::vector<stk::mesh::Selector>& selectors, const stk::mesh::impl::LocalIdMapper& localIds)
{
    Teuchos::ParameterList params = getGraphBasedParameters(balanceSettings, numSubdomainsToCreate);

    std::ostringstream os;
    os << "Using Zoltan2 version: " << Zoltan2::Zoltan2_Version();
    logMessage(stkMeshBulkData.parallel(), os.str());
    logMessage(stkMeshBulkData.parallel(), "Filling in graph data");


    Zoltan2ParallelGraph zoltan2Graph;
    createZoltanParallelGraph(balanceSettings, stkMeshBulkData, selectors, localIds, zoltan2Graph);

    std::vector<double> copyOrigWeights = zoltan2Graph.get_vertex_weights();
    std::vector<int> all_local_ids(copyOrigWeights.size());
    for(size_t i=0;i<all_local_ids.size();++i)
    {
        all_local_ids[i]=i;
    }

    decomp.clear();
    decomp.resize(copyOrigWeights.size(), std::make_pair(stk::mesh::Entity(), stkMeshBulkData.parallel_rank()));

    #if defined(WRITE_OUT_DEBUGGING_INFO) || defined(WRITE_OUT_DECOMP_METRICS)
    static int step = 0;
    stk::mesh::EntityVector entities = save_for_debugging_local_ids(stkMeshBulkData);
    #endif

    if (balanceSettings.isMultiCriteriaRebalance())
    {
        stk::mesh::Selector selectUnion = stk::mesh::selectUnion(selectors);
        get_multicriteria_parmetis_decomp(stkMeshBulkData, balanceSettings, zoltan2Graph, params, selectUnion, decomp, localIds);
    }
    else
    {
        for(size_t i=0;i<selectors.size();++i)
        {
            zoltan2Graph.set_vertex_weights(copyOrigWeights);

            std::vector<int> localIdsOfNotSelectedEntities = getLocalIdsOfEntitiesNotSelected(stkMeshBulkData, selectors[i], localIds);

            size_t numItemsToRebalance = copyOrigWeights.size() - localIdsOfNotSelectedEntities.size();
            size_t numVerticesGlobal = 0;
            stk::all_reduce_sum(stkMeshBulkData.parallel(), &numItemsToRebalance, &numVerticesGlobal, 1);

            if(numVerticesGlobal>0)
            {
                for(size_t j=0;j<localIdsOfNotSelectedEntities.size();++j)
                {
                    zoltan2Graph.set_vertex_weight(localIdsOfNotSelectedEntities[j], 0.0);
                }

                get_multicriteria_parmetis_decomp(stkMeshBulkData, balanceSettings, zoltan2Graph, params, selectors[i], decomp, localIds);
            }
        }
    }

    logMessage(stkMeshBulkData.parallel(), "Finished decomposition solve");

    #if defined(WRITE_OUT_DECOMP_METRICS)
    std::vector<double> weights_per_proc(numSubdomainsToCreate,0);
    double max = 0;
    double sum = 0;
    for(size_t i=0; i<zoltan2Graph.mVertexIds.size();++i)
    {
        sum += copyOrigWeights[i];
        max = std::max(max, copyOrigWeights[i]);
        weights_per_proc[decomp[i]] += copyOrigWeights[i];
    }

    double global_sum = 0;
    stk::all_reduce_sum(stkMeshBulkData.parallel(), &sum, &global_sum, 1);
    double global_max = 0;
    stk::all_reduce_max(stkMeshBulkData.parallel(), &max, &global_max, 1);
    double max_weight_any_proc = 0;
    stk::all_reduce_max(stkMeshBulkData.parallel(), &sum, &max_weight_any_proc, 1);

    std::vector<double> weights_per_proc_sum(weights_per_proc.size(),0);
    stk::all_reduce_sum(stkMeshBulkData.parallel(), weights_per_proc.data(), weights_per_proc_sum.data(), weights_per_proc.size());

    {
        std::ostringstream filename;
        filename << "balance_info_" << stkMeshBulkData.parallel_rank() << ".txt";
        std::ofstream out(filename.str().c_str(), std::ofstream::app);
        std::ostringstream os;
        os << "=========================== For Step " << step << " ====================================" << std::endl;
        os << "Max element weight: " << global_max << " and average weight per proc: " << global_sum/stkMeshBulkData.parallel_size() << " and max total weight any proc = " << max_weight_any_proc << std::endl;
        os << "Max element weight/Ave weight per proc: " << global_max/std::max(0.000001, global_sum/stkMeshBulkData.parallel_size()) << std::endl;
        os << "Weight before this proc: " << sum << " and weight after decomp: " << weights_per_proc_sum[stkMeshBulkData.parallel_rank()] << std::endl;
        os << "Weight/average before: " << sum/(std::max(0.000001, global_sum/stkMeshBulkData.parallel_size())) << " and after: " << weights_per_proc_sum[stkMeshBulkData.parallel_rank()]/(std::max(0.000001, global_sum/stkMeshBulkData.parallel_size())) << std::endl;
        out << os.str();
        out.close();
    }
    #endif

    #if defined(WRITE_OUT_DEBUGGING_INFO) || defined(WRITE_OUT_DECOMP_METRICS)
    step++;
    #endif
}

void calculateGeometricOrGraphBasedDecomp(const BalanceSettings& balanceSettings, const int numSubdomainsToCreate, stk::mesh::EntityProcVec &decomp, stk::mesh::BulkData& stkMeshBulkData, const std::vector<stk::mesh::Selector>& selectors)
{
    ThrowRequireWithSierraHelpMsg(numSubdomainsToCreate > 0);
    ThrowRequireWithSierraHelpMsg(is_geometric_method(balanceSettings.getDecompMethod()) || balanceSettings.getDecompMethod()=="parmetis" || balanceSettings.getDecompMethod()=="zoltan");

    if(is_geometric_method(balanceSettings.getDecompMethod()))
    {
        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);
        fill_decomp_using_geometric_method(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors, localIds);
    }
    else if (balanceSettings.getDecompMethod()=="parmetis" || balanceSettings.getDecompMethod()=="zoltan")
    {
        stk::mesh::Ghosting * customAura = internal::create_custom_ghosting(stkMeshBulkData, balanceSettings);
        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);
        internal::fill_connectivity_count_field(stkMeshBulkData, balanceSettings);
        fill_decomp_using_parmetis(balanceSettings, numSubdomainsToCreate, decomp, stkMeshBulkData, selectors, localIds);
        internal::destroy_custom_ghosting(stkMeshBulkData, customAura);
    }
}

bool compareEntityEqualityOnly(const std::pair<stk::mesh::Entity,int> &a, const std::pair<stk::mesh::Entity,int> &b)
{
    return (a.first == b.first);
}

void add_if_owned(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity, int newOwningProc, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    if(stkMeshBulkData.bucket(entity).owned())
        entityProcPairs.emplace_back(entity, newOwningProc);
}

void add_connected_entities_of_rank(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element, int newOwningProc, stk::mesh::EntityRank rank, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    unsigned numEntities = stkMeshBulkData.num_connectivity(element, rank);
    const stk::mesh::Entity *entities = stkMeshBulkData.begin(element, rank);
    for(unsigned int i = 0; i < numEntities; i++)
    {
        add_if_owned(stkMeshBulkData, entities[i], newOwningProc, entityProcPairs);
    }
}

void fillEntityProcPairsForEntityAndItsClosure(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity elementToMove, int newOwningProc,
        std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    entityProcPairs.emplace_back(elementToMove, newOwningProc);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::FACE_RANK, entityProcPairs);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::EDGE_RANK, entityProcPairs);
    add_connected_entities_of_rank(stkMeshBulkData, elementToMove, newOwningProc, stk::topology::NODE_RANK, entityProcPairs);
}

void performModifications(stk::mesh::BulkData& stkMeshBulkData, std::vector<std::pair<stk::mesh::Entity, int> > &entityProcPairs)
{
    std::sort(entityProcPairs.begin(), entityProcPairs.end());
    std::vector<std::pair<stk::mesh::Entity, int> >::iterator iter = std::unique(entityProcPairs.begin(), entityProcPairs.end(), compareEntityEqualityOnly);
    entityProcPairs.erase(iter, entityProcPairs.end());

    stkMeshBulkData.change_entity_owner(entityProcPairs);
}

void rebalance(DecompositionChangeList &changeList)
{
    stk::mesh::EntityProcVec entityProcPairs;
    changeList.pImpl->get_decomposition_with_full_closure(entityProcPairs);
    performModifications(changeList.get_bulk(), entityProcPairs);
}

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec entityProcPairs;
    for(const stk::mesh::EntityProc& entity_proc : decomposition)
    {
        if(entity_proc.second != stkMeshBulkData.parallel_rank())
            fillEntityProcPairsForEntityAndItsClosure(stkMeshBulkData, entity_proc.first, entity_proc.second, entityProcPairs);
    }
    performModifications(stkMeshBulkData, entityProcPairs);
}

stk::mesh::EntityProcVec get_mapped_decomp(const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec mapped_decomp(decomposition.size());
    for(size_t i=0;i<decomposition.size();++i)
    {
        mapped_decomp[i].first = decomposition[i].first;
        mapped_decomp[i].second = mappings[decomposition[i].second];
    }
    return mapped_decomp;
}

void rebalance(stk::mesh::BulkData& stkMeshBulkData, const std::vector<unsigned>& mappings, const stk::mesh::EntityProcVec& decomposition)
{
    stk::mesh::EntityProcVec mapped_decomp = get_mapped_decomp(mappings, decomposition);
    rebalance(stkMeshBulkData, mapped_decomp);
}

void print_rebalance_metrics(const size_t num_global_entity_migrations, const size_t max_global_entity_migrations, stk::mesh::BulkData & stkMeshBulkData)
{
    std::vector<size_t> meshCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, meshCounts);
    double fraction_of_mesh_moved = static_cast<double>(num_global_entity_migrations)/meshCounts[stk::topology::ELEM_RANK];
    double avg_global_entity_migrations = static_cast<double>(num_global_entity_migrations)/stkMeshBulkData.parallel_size();

    std::ostringstream oss;
    oss << "Percentage of total mesh moved = " << fraction_of_mesh_moved*100.0 << "%";
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(),oss.str());

    oss.str("");
    oss << "Max/Avg global entity migrations = " << max_global_entity_migrations/avg_global_entity_migrations;
    stk::balance::internal::logMessage(stkMeshBulkData.parallel(),oss.str());
}

} //internal

}
}

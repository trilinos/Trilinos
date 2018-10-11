#include "DetectAndFixMechanisms.hpp"

#include "Zoltan2ParallelGraph.hpp"
#include "balanceTypes.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "stk_balance/internal/privateDeclarations.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_util/parallel/CommSparse.hpp"

namespace stk { namespace balance { namespace internal {

void translate_data_for_dsd_detection(const Zoltan2ParallelGraph &zoltan2Graph, const stk::mesh::BulkData& bulk, const stk::mesh::impl::LocalIdMapper& localIds,
                                      std::vector<BalanceLocalNumber> &tmpOffsets, std::vector<BalanceGlobalNumber>& tmpAdj);

int extract_connected_lists( const std::vector<BalanceGlobalNumber>& adj,
                             const std::vector<BalanceLocalNumber>& offsets,
                             std::vector<BalanceGlobalNumber>& list,
                             std::vector<BalanceLocalNumber>& list_ptr );

std::vector<stk::mesh::EntityVector> get_elements_per_component(const stk::mesh::BulkData& bulk, const std::vector<BalanceLocalNumber> &list_ptr, const std::vector<BalanceGlobalNumber>& list);

std::vector<int> get_components_to_move(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::EntityVector>& elementsPerComponent);

bool detectAndFixMechanisms(const stk::balance::BalanceSettings& graphSettings, stk::mesh::BulkData &bulk)
{
    stk::mesh::Ghosting * customAura = internal::create_custom_ghosting(bulk, graphSettings);
    stk::mesh::impl::LocalIdMapper localIds(bulk, stk::topology::ELEM_RANK);

    Zoltan2ParallelGraph zoltan2Graph;
    fill_zoltan2_parallel_graph(bulk, graphSettings, localIds, zoltan2Graph);

    std::vector<BalanceLocalNumber> tmpOffsets;
    std::vector<BalanceGlobalNumber> tmpAdj;
    translate_data_for_dsd_detection(zoltan2Graph, bulk, localIds, tmpOffsets, tmpAdj);

    std::vector<BalanceGlobalNumber> list;
    std::vector<BalanceLocalNumber> list_ptr;
    size_t num_dsd = extract_connected_lists(tmpAdj, tmpOffsets, list, list_ptr);

    std::vector<stk::mesh::EntityVector> elementsPerComponent;
    std::vector<int> componentsToMove;
    if(num_dsd > 1)
    {
        ThrowRequireWithSierraHelpMsg(num_dsd == list_ptr.size()-1);
        elementsPerComponent = get_elements_per_component(bulk, list_ptr, list);
        componentsToMove = get_components_to_move(bulk, elementsPerComponent);
    }

    bool hasMechanisms = elementsPerComponent.size()>1;
    bool globallyHaveMechanisms = stk::is_true_on_any_proc(bulk.parallel(), hasMechanisms);

    if(globallyHaveMechanisms) {
        move_components(zoltan2Graph, localIds, bulk, elementsPerComponent, componentsToMove);
    }

    internal::destroy_custom_ghosting(bulk, customAura);

    return globallyHaveMechanisms;
}


void fill_zoltan2_parallel_graph(stk::mesh::BulkData& bulk, const stk::balance::BalanceSettings& graphSettings, const stk::mesh::impl::LocalIdMapper& localIds, Zoltan2ParallelGraph &zoltan2Graph)
{
    stk::mesh::Selector sel = bulk.mesh_meta_data().locally_owned_part();
    zoltan2Graph.setMechanismCheckFlag();
    stk::balance::internal::createZoltanParallelGraph(graphSettings, bulk, std::vector<stk::mesh::Selector>{sel}, localIds, zoltan2Graph);
}

size_t one_based(size_t input)
{
    return input+1;
}

void translate_data_for_dsd_detection(const Zoltan2ParallelGraph &zoltan2Graph, const stk::mesh::BulkData& bulk, const stk::mesh::impl::LocalIdMapper& localIds,
                                      std::vector<BalanceLocalNumber> &tmpOffsets, std::vector<BalanceGlobalNumber>& tmpAdj)
{
    const std::vector<BalanceLocalNumber> &offsets = zoltan2Graph.get_offsets();
    const std::vector<BalanceGlobalNumber> &adj = zoltan2Graph.get_adjacency();

    tmpOffsets.resize(offsets.size());

    tmpOffsets[0] = one_based(0);
    int numOffsets = 0;

    for(size_t i=0;i<offsets.size()-1;++i)
    {
        // start with self connectivity
        numOffsets++;
        tmpAdj.push_back(one_based(i));

        size_t start = offsets[i];
        size_t end = offsets[i+1];
        for(size_t j=start;j<end;++j)
        {
            stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, adj[j]);
            if(bulk.is_valid(element) && bulk.bucket(element).owned())
            {
                numOffsets++;
                tmpAdj.push_back(one_based(localIds.entity_to_local(element)));
            }
        }
        tmpOffsets[i+1] = one_based(numOffsets);
    }
}

std::set<stk::mesh::Entity> getNodesOfComponent(const stk::mesh::BulkData& bulk, const stk::mesh::EntityVector& elements)
{
    std::set<stk::mesh::Entity> nodesOfComponent;
    for(size_t i=0;i<elements.size();++i)
    {
        stk::mesh::Entity element = elements[i];
        nodesOfComponent.insert(bulk.begin_nodes(element), bulk.end_nodes(element));
    }
    return nodesOfComponent;
}


std::vector<stk::mesh::EntityVector> get_elements_per_component(const stk::mesh::BulkData& bulk, const std::vector<BalanceLocalNumber> &list_ptr, const std::vector<BalanceGlobalNumber>& list)
{
    size_t num_dsd = list_ptr.size()-1;
    ThrowRequireWithSierraHelpMsg(num_dsd>1);

    std::vector<stk::mesh::EntityVector> elementsPerComponent(num_dsd);

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK), elements);

    for(size_t i=0;i<list_ptr.size()-1;++i)
    {
        size_t start = list_ptr[i];
        size_t end = list_ptr[i+1];
        for(size_t j=start;j<end;++j)
        {
            stk::mesh::Entity element = elements[list[j]-1];
            elementsPerComponent[i].push_back(element);
        }
    }

    return elementsPerComponent;
}

size_t get_max_component_index(const std::vector<stk::mesh::EntityVector> &elementsPerComponent)
{
    size_t max_component_index = 0;
    size_t max_elements = 0;
    for(size_t i=0;i<elementsPerComponent.size();++i)
    {
        if(elementsPerComponent[i].size() > max_elements)
        {
            max_component_index = i;
            max_elements = elementsPerComponent[i].size();
        }
    }
    return max_component_index;
}

stk::mesh::EntityVector get_shared_nodes_between_components(const std::set<stk::mesh::Entity>& nodesOfHomeComponent, const std::set<stk::mesh::Entity>& nodesOfThisComponent)
{
    std::vector<stk::mesh::Entity> result(std::max(nodesOfHomeComponent.size(), nodesOfThisComponent.size()));
    std::vector<stk::mesh::Entity>::iterator iter=std::set_intersection(nodesOfHomeComponent.begin(), nodesOfHomeComponent.end(),
                                                                        nodesOfThisComponent.begin(), nodesOfThisComponent.end(),
                                                                        result.begin());

    result.resize(iter-result.begin());
    return result;
}

std::vector<int> get_components_to_move(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::EntityVector>& elementsPerComponent)
{
    size_t maxComponent = get_max_component_index(elementsPerComponent);
    std::vector<int> componentsToMove;

    std::set<stk::mesh::Entity> nodesOfHomeComponent = getNodesOfComponent(bulk, elementsPerComponent[maxComponent]);

    std::ostringstream os;
    os << "Proc " << bulk.parallel_rank() << std::endl;

    for(size_t i=0;i<elementsPerComponent.size();++i)
    {
        if(i!=maxComponent)
        {
            std::set<stk::mesh::Entity> nodesOfThisComponent = getNodesOfComponent(bulk, elementsPerComponent[i]);
            stk::mesh::EntityVector sharedNodesBetweenComponents = get_shared_nodes_between_components(nodesOfHomeComponent, nodesOfThisComponent);

            if(sharedNodesBetweenComponents.size()>0)
            {
                os << "found a mechanism between components on nodes: ";
                for(size_t j=0;j<sharedNodesBetweenComponents.size();++j)
                    os << bulk.identifier(sharedNodesBetweenComponents[j]) << " ";
                os << std::endl;

                componentsToMove.push_back(i);
            }
            else
            {
                nodesOfHomeComponent.insert(nodesOfThisComponent.begin(), nodesOfThisComponent.end());
            }
        }
    }

    if(componentsToMove.size()>=1)
        std::cerr << os.str();

    return componentsToMove;
}



std::vector<stk::mesh::EntityId> get_elements_on_the_move(const stk::mesh::BulkData &bulk, const std::set<stk::mesh::EntityProc>& entityProcs)
{
    std::vector<stk::mesh::EntityId> elementsOnTheMove;
    stk::CommSparse comm(bulk.parallel());

    std::set<stk::mesh::EntityProc>::iterator iter;
    for(int iphase=0;iphase<2;++iphase)
    {
        for(iter=entityProcs.begin();iter!=entityProcs.end();++iter)
        {
            stk::mesh::EntityId id = bulk.identifier(iter->first);
            int proc = iter->second;
            comm.send_buffer(proc).pack<stk::mesh::EntityId>(id);
        }

        if(iphase==0)
            comm.allocate_buffers();
        else
            comm.communicate();
    }

    for(int proc=0;proc<bulk.parallel_size();++proc)
    {
        stk::CommBuffer & buf = comm.recv_buffer(proc);
        while(buf.remaining())
        {
            stk::mesh::EntityId id;
            buf.unpack(id);
            elementsOnTheMove.push_back(id);
        }
    }
    return elementsOnTheMove;
}

stk::mesh::EntityProcVec get_element_proc_movement(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::EntityId>& otherElementsOnTheMove,
                                                   const std::vector<stk::mesh::EntityVector>& elementsToMove,
                                                   const stk::mesh::impl::LocalIdMapper localIds,
                                                   const std::vector<BalanceLocalNumber> &offsets, const std::vector<BalanceGlobalNumber> &adj)
{
    stk::mesh::EntityProcVec element_migration;

    for(size_t i=0;i<elementsToMove.size();++i)
    {
        int destProcForComp = -1;
        for(size_t j=0;j<elementsToMove[i].size();++j)
        {
            stk::mesh::Entity element = elementsToMove[i][j];
            unsigned local_id = localIds.entity_to_local(element);
            size_t start = offsets[local_id];
            size_t end   = offsets[local_id+1];
            for(size_t k=start;k<end;++k)
            {
                stk::mesh::Entity otherElement = bulk.get_entity(stk::topology::ELEM_RANK, adj[k]);
                if(bulk.is_valid(otherElement) && !bulk.bucket(otherElement).owned())
                {
                    stk::mesh::EntityId id = bulk.identifier(otherElement);
                    std::vector<stk::mesh::EntityId>::const_iterator iter = std::find(otherElementsOnTheMove.begin(),
                                                                                      otherElementsOnTheMove.end(),
                                                                                      id);
                    if(iter == otherElementsOnTheMove.end())
                    {
                        std::vector<int> procs;
                        bulk.comm_procs(bulk.entity_key(otherElement), procs);
                        ThrowRequireWithSierraHelpMsg(procs.size()>0);
                        destProcForComp = procs[0];
                        break;
                    }
                }
            }

            if(destProcForComp != -1)
                break;
        }

        if(destProcForComp == -1)
        {
            if(bulk.parallel_rank() == 0)
                std::cerr << "Global mesh has mechanism!" << std::endl;
        }
        else
        {
            for(size_t j=0;j<elementsToMove[i].size();++j)
            {
                stk::mesh::Entity element = elementsToMove[i][j];
                if(bulk.is_valid(element) && bulk.bucket(element).owned())
                    element_migration.push_back(stk::mesh::EntityProc(element, destProcForComp));
            }
        }
    }

    return element_migration;
}

std::set<stk::mesh::EntityProc> get_element_proc_info(const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::EntityVector>& elementsToMove,
                                                      const stk::mesh::impl::LocalIdMapper localIds,
                                                      const std::vector<BalanceLocalNumber> &offsets, const std::vector<BalanceGlobalNumber> &adj)
{
    std::set<stk::mesh::EntityProc> entityProcs;

    for(size_t i=0;i<elementsToMove.size();++i)
    {
        for(size_t j=0;j<elementsToMove[i].size();++j)
        {
            stk::mesh::Entity element = elementsToMove[i][j];
            unsigned local_id = localIds.entity_to_local(element);
            size_t start = offsets[local_id];
            size_t end   = offsets[local_id+1];
            for(size_t k=start;k<end;++k)
            {
                stk::mesh::Entity otherElement = bulk.get_entity(stk::topology::ELEM_RANK, adj[k]);
                if(bulk.is_valid(otherElement) && !bulk.bucket(otherElement).owned())
                {
                    std::vector<int> procs;
                    bulk.comm_procs(bulk.entity_key(otherElement), procs);
                    for(size_t m=0;m<procs.size();++m)
                    {
                        entityProcs.insert(stk::mesh::EntityProc(element, procs[m]));
                    }
                }
            }
        }
    }

    return entityProcs;
}

void move_components(const Zoltan2ParallelGraph &zoltan2Graph,
                     const stk::mesh::impl::LocalIdMapper &localIds,
                     stk::mesh::BulkData& bulk,
                     const std::vector<stk::mesh::EntityVector>& elementsToMove,
                     const std::vector<int>& componentsToMove)
{


    const std::vector<BalanceLocalNumber> &offsets = zoltan2Graph.get_offsets();
    const std::vector<BalanceGlobalNumber> &adj = zoltan2Graph.get_adjacency();

    std::vector<stk::mesh::EntityVector> localCopyElementsToMove;
    localCopyElementsToMove.reserve(elementsToMove.size());

    for(size_t i=0;i<componentsToMove.size();++i)
        localCopyElementsToMove.push_back(elementsToMove[componentsToMove[i]]);

    std::set<stk::mesh::EntityProc> entityProcs = get_element_proc_info(bulk, localCopyElementsToMove, localIds, offsets, adj);

    std::vector<stk::mesh::EntityId> elementsFromOtherProcsOnTheMove = get_elements_on_the_move(bulk, entityProcs);
    std::sort(elementsFromOtherProcsOnTheMove.begin(), elementsFromOtherProcsOnTheMove.end());

    stk::mesh::EntityProcVec element_migration = get_element_proc_movement(bulk, elementsFromOtherProcsOnTheMove, localCopyElementsToMove, localIds, offsets, adj);
    DecompositionChangeList changeList(bulk, element_migration);

    const size_t num_global_entity_migrations = changeList.get_num_global_entity_migrations();
    if(num_global_entity_migrations>0)
        internal::rebalance(changeList);
}

int extract_connected_lists( const std::vector<BalanceGlobalNumber>& adj,
                             const std::vector<BalanceLocalNumber>& offsets,
                             std::vector<BalanceGlobalNumber>& list,
                             std::vector<BalanceLocalNumber>& list_ptr )
{
    int components = 1;
    size_t nrow = offsets.size()-1;
    if(nrow != 0)
    {
        int root = 1;
        std::vector<int> mask(nrow, 1);
        size_t nordered = 1;
        list_ptr.push_back(nordered - 1);
        list.push_back(root);
        mask[root - 1] = 0;

        int ni = 1;
        int nf = 1;

        int nni = 0, nnf = 0;

        while(nordered < nrow)
        {
            if(nf == ni - 1)
            {
                ++components;
                list_ptr.push_back(nordered);
                for(size_t i = 0; i < nrow; i++)
                {
                    if(mask[i] == 1)
                    {
                        ++nordered;
                        list.push_back(i + 1);
                        mask[i] = 0;
                        nf = ni;
                        break;
                    }
                }
            }
            nni = nf + 1;
            nnf = nni - 1;
            for(int i = ni; i <= nf; i++)
            {
                int ki = offsets[list[i - 1] - 1] - 1;
                int kf = offsets[list[i - 1]] - 2;
                for(int k = ki; k <= kf; k++)
                {
                    int j = adj[k];
                    if(mask[j - 1] == 1)
                    {
                        ++nnf;
                        ++nordered;
                        mask[j - 1] = 0;
                        list.push_back(j);
                    }
                }
            }
            ni = nni;
            nf = nnf;
        }
        list_ptr.push_back(nordered);
    }
    else
        components = 0;

    return components;
}


}}} // end of namespaces

#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "../../stk_util/stk_util/parallel/DistributedIndex.hpp"
#include "../../stk_mesh/stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "../../stk_mesh/stk_mesh/base/GetEntities.hpp"
#include "../../stk_util/stk_util/parallel/ParallelReduce.hpp"
#include "../../stk_mesh/stk_mesh/base/FEMHelpers.hpp"

namespace stk { namespace mesh {
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

void print_and_throw_if_entities_are_not_unique(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs)
{
    bool allOk = badKeyProcs.empty();
    if(!allOk)
        std::cerr << get_non_unique_key_messages(bulkData, badKeyProcs);

    bool globally_ok = stk::is_true_on_all_procs(bulkData.parallel(), allOk);
    ThrowRequireMsg(globally_ok, "Program error. Please contact sierra-help@sandia.gov for support.");
}

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData)
{
    stk::parallel::DistributedIndex distributedIndex( bulkData.parallel(), stk::mesh::impl::convert_entity_keys_to_spans(bulkData.mesh_meta_data()));
    add_keys_to_distributed_index(bulkData, distributedIndex);
    stk::parallel::DistributedIndex::KeyTypeVector localKeys = get_all_local_keys(bulkData);
    return get_non_unique_keys(bulkData, distributedIndex, localKeys);
}

}}

class MeshChecker : public stk::unit_test_util::MeshFixture
{
public:
};

void test_non_unique_key_results(const stk::mesh::BulkData& bulkData, const std::vector<std::vector<stk::mesh::EntityKeyProc>>& gold_values)
{
    std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(bulkData);
    EXPECT_EQ(gold_values[bulkData.parallel_rank()].size(), badKeyProcs.size());
    EXPECT_EQ(gold_values[bulkData.parallel_rank()], badKeyProcs);
    EXPECT_THROW(stk::mesh::print_and_throw_if_entities_are_not_unique(bulkData, badKeyProcs), std::logic_error);
}

void test_keys_are_unique(const stk::mesh::BulkData& bulkData)
{
    std::vector<stk::mesh::EntityKeyProc> badKeyProcs = stk::mesh::get_non_unique_key_procs(bulkData);
    EXPECT_TRUE(badKeyProcs.empty());
    EXPECT_NO_THROW(stk::mesh::print_and_throw_if_entities_are_not_unique(bulkData, badKeyProcs));
}

TEST_F(MeshChecker, check_throw_on_non_unique_ids)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0}},
        };

        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        get_bulk().modification_begin();
        get_bulk().declare_entity(stk::topology::NODE_RANK, 1);
        get_bulk().modification_end();

        test_non_unique_key_results(get_bulk(), gold_values);
    }
}

TEST_F(MeshChecker, check_no_throw_on_shared_nodes)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        get_bulk().modification_begin();
        stk::mesh::Entity node1 = get_bulk().declare_entity(stk::topology::NODE_RANK, 1);
        get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
        get_bulk().modification_end();

        test_keys_are_unique(get_bulk());
    }
}

TEST_F(MeshChecker, check_throw_on_partial_shared_nodes)
{
    if(stk::parallel_machine_size(get_comm())==3)
    {
        std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
                 { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
        };

        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        get_bulk().modification_begin();
        stk::mesh::Entity node1 = get_bulk().declare_entity(stk::topology::NODE_RANK, 1);
        if(get_bulk().parallel_rank()<2)
        {
            get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
        }
        get_bulk().modification_end();

        test_non_unique_key_results(get_bulk(), gold_values);
    }
}

TEST_F(MeshChecker, check_throw_on_partial_shared_nodes_on_four_procs)
{
    if(stk::parallel_machine_size(get_comm())==4)
    {
        std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2},
                 { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 3}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 2},
                 { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 3}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
                 { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
                {{ stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 0},
                 { stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), 1}},
        };

        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        get_bulk().modification_begin();
        stk::mesh::Entity node1 = get_bulk().declare_entity(stk::topology::NODE_RANK, 1);
        if(get_bulk().parallel_rank()<2)
        {
            get_bulk().add_node_sharing(node1, 1-get_bulk().parallel_rank());
        }
        else
        {
            const int other_proc = get_bulk().parallel_rank()==2 ? 3 : 2;
            get_bulk().add_node_sharing(node1, other_proc);
        }
        get_bulk().modification_end();

        test_non_unique_key_results(get_bulk(), gold_values);
    }
}

TEST_F(MeshChecker, check_throw_on_non_unique_ids_for_elements)
{
    if(stk::parallel_machine_size(get_comm())==2)
    {
        std::vector<std::vector<stk::mesh::EntityKeyProc>> gold_values = {
                {{ stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 1}},
                {{ stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1), 0}},
        };

        setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
        get_bulk().modification_begin();
        stk::mesh::EntityId nodeId = static_cast<stk::mesh::EntityId>(get_bulk().parallel_rank()+1u);
        stk::mesh::declare_element(get_bulk(),get_meta().get_topology_root_part(stk::topology::PARTICLE),1u,stk::mesh::EntityIdVector{nodeId});
        get_bulk().modification_end();

        test_non_unique_key_results(get_bulk(), gold_values);
    }
}


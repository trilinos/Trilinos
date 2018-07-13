#include "MtoNRebalancer.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/entityDataToField.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>

namespace stk {
namespace balance {
namespace internal {

MtoNRebalancer::MtoNRebalancer(stk::mesh::BulkData &bulkData, stk::mesh::Field<double> &targetField,
                               const stk::balance::BalanceSettings &graphSettings, int num_target_procs)
: mBulkData(bulkData), targetDecompField(targetField), mNumTargetProcs(num_target_procs), decomp()
{
    decomp = stk::balance::internal::get_element_decomp(mNumTargetProcs, mBulkData, graphSettings);
}

MtoNRebalancer::~MtoNRebalancer() {}

void MtoNRebalancer::move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc)
{
    store_off_target_proc_on_elements_before_moving_subdomains();
    stk::balance::internal::rebalance(mBulkData, target_proc_to_starting_proc, decomp);
    move_entities_into_mapped_subdomain_parts(target_proc_to_starting_proc);
}

void MtoNRebalancer::move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings)
{
    declare_all_subdomain_parts();
    mBulkData.modification_begin();
    change_parts_on_entities_on_all_subdomains(mappings);
    mBulkData.modification_end();
}

void MtoNRebalancer::change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping)
{
    for(int i = 0; i < mNumTargetProcs; i++)
    {
        if(subdomain_proc_mapping[i] == static_cast<unsigned>(mBulkData.parallel_rank()))
        {
            std::vector<stk::mesh::Entity> entities = get_entities_for_subdomain(i);
            move_entities_into_subdomain_part(i, entities);
        }
    }
}

stk::mesh::BulkData& MtoNRebalancer::get_bulk()
{
    return mBulkData;
}

void MtoNRebalancer::move_entities_into_subdomain_part(size_t i, const stk::mesh::EntityVector &entities)
{
    stk::mesh::PartVector partVector = get_parts_to_add_for_subdomain(i);
    for(size_t j = 0; j < entities.size(); j++)
        mBulkData.change_entity_parts(entities[j], partVector);
}

stk::mesh::PartVector MtoNRebalancer::get_parts_to_add_for_subdomain(size_t subdomain_num)
{
    stk::mesh::Part* part = get_subdomain_part(subdomain_num);
    return {part};
}

stk::mesh::Part* MtoNRebalancer::get_subdomain_part(size_t subdomain_num)
{
    std::string partNameForSubdomain = getSubdomainPartName(subdomain_num);
    return get_meta().get_part(partNameForSubdomain);
}

std::string MtoNRebalancer::getSubdomainPartName(int subdomainId)
{
    std::ostringstream out;
    out << "subdomain_" << subdomainId;
    return out.str();
}

std::vector<stk::mesh::Entity> MtoNRebalancer::get_entities_for_subdomain(size_t subdomain_num)
{
    const stk::mesh::BucketVector &buckets = get_bulk().buckets(stk::topology::ELEMENT_RANK);
    return get_entitites_for_subdomain_using_field_from_buckets(subdomain_num, buckets);
}

stk::mesh::EntityVector MtoNRebalancer::get_entitites_for_subdomain_using_field_from_buckets(size_t subdomain_num, const stk::mesh::BucketVector& buckets)
{
    std::vector<stk::mesh::Entity> entities;
    for(size_t j = 0; j < buckets.size(); j++)
        add_owned_entities_from_bucket_using_target_decomp_field(*buckets[j], subdomain_num, entities);
    return entities;
}

void MtoNRebalancer::add_owned_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities)
{
    if(bucket.owned())
        add_entities_from_bucket_using_target_decomp_field(bucket, subdomain_num, entities);
}

void MtoNRebalancer::add_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities)
{
    double *bucketSubdomainData = static_cast<double*>(stk::mesh::field_data(targetDecompField, bucket));
    for(size_t k = 0; k < bucket.size(); k++)
    {
        if(bucketSubdomainData[k] == static_cast<double>(subdomain_num))
            entities.push_back(bucket[k]);
    }
}

stk::mesh::MetaData& MtoNRebalancer::get_meta()
{
    return mBulkData.mesh_meta_data();
}

void MtoNRebalancer::declare_all_subdomain_parts()
{
    for(int i=0;i<mNumTargetProcs;++i)
    {
        std::string partNameForSubdomain = getSubdomainPartName(i);
        get_meta().declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK);
    }
}

void MtoNRebalancer::store_off_target_proc_on_elements_before_moving_subdomains()
{
    stk::balance::internal::put_entity_data_to_field(decomp, &targetDecompField);
}

bool MtoNRebalancer::does_this_proc_own_subdomain(unsigned subdomainOwner)
{
    return (subdomainOwner == static_cast<unsigned>(get_bulk().parallel_rank()));
}

void MtoNRebalancer::create_subdomain_and_write(const std::string &filename,
                                                unsigned subdomain,
                                                int global_num_nodes,
                                                int global_num_elems,
                                                const stk::io::EntitySharingInfo &nodeSharingInfo,
                                                int numSteps, double timeStep)
{
    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
    stk::tools::copy_mesh(get_bulk(), *get_meta().get_part(getSubdomainPartName(subdomain)), newBulkData);
    stk::io::write_file_for_subdomain(filename, subdomain, mNumTargetProcs, global_num_nodes, global_num_elems, newBulkData, nodeSharingInfo, numSteps, timeStep);
}

stk::io::EntitySharingInfo MtoNRebalancer::get_node_sharing_info(unsigned subdomain)
{
    stk::mesh::EntityVector sharedNodes;
    std::vector<int> procsForSharedNodes;
    fill_shared_node_info_for_this_subdomain(subdomain, sharedNodes, procsForSharedNodes);

    stk::io::EntitySharingInfo nodeSharingInfo;
    for(size_t nodeIndex = 0; nodeIndex < sharedNodes.size(); nodeIndex++)
    {
        stk::mesh::EntityId nodeId = get_bulk().identifier(sharedNodes[nodeIndex]);
        nodeSharingInfo.push_back({nodeId, procsForSharedNodes[nodeIndex]});
    }
    return nodeSharingInfo;
}

stk::mesh::EntityVector MtoNRebalancer::get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index)
{
    stk::mesh::Selector selected_nodes = *get_subdomain_part(this_subdomain_index) & *get_subdomain_part(other_subdomain_index);
    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(selected_nodes, get_bulk().buckets(stk::topology::NODE_RANK), nodes);
    return nodes;
}

void MtoNRebalancer::fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num)
{
    stk::mesh::EntityVector nodes = get_nodes_shared_between_subdomains(this_subdomain_num, other_subdomain_num);
    shared_nodes.insert(shared_nodes.end(), nodes.begin(), nodes.end());
    procs_for_shared_nodes.resize(shared_nodes.size(), other_subdomain_num);
}

void MtoNRebalancer::fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes)
{
    for(int other_subdomain_num=0;other_subdomain_num<mNumTargetProcs;++other_subdomain_num)
    {
        if(static_cast<int>(this_subdomain_num) != other_subdomain_num)
            fill_shared_node_proc_info(shared_nodes, procs_for_shared_nodes, this_subdomain_num, other_subdomain_num);
    }
}

}}}

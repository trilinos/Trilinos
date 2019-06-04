#include "MtoNRebalancer.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/entityDataToField.hpp>
#include <stk_balance/internal/MxNutils.hpp>

namespace stk {
namespace balance {
namespace internal {

MtoNRebalancer::MtoNRebalancer(stk::mesh::BulkData &bulkData, stk::mesh::Field<double> &targetField,
                               const stk::balance::BalanceSettings &graphSettings, int numTargetProcs)
: mBulkData(bulkData), subdomainCreator(bulkData, numTargetProcs), targetDecompField(targetField), decomp()
{
    decomp = stk::balance::internal::get_element_decomp(subdomainCreator.get_num_subdomains(), mBulkData, graphSettings);
}

MtoNRebalancer::~MtoNRebalancer() {}

void MtoNRebalancer::move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc)
{
    store_off_target_proc_on_elements_before_moving_subdomains();
    stk::balance::internal::rebalance(mBulkData, target_proc_to_starting_proc, decomp);
    move_entities_into_mapped_subdomain_parts(target_proc_to_starting_proc);
}

stk::io::EntitySharingInfo MtoNRebalancer::get_node_sharing_info(unsigned subdomain)
{
    return subdomainCreator.get_node_sharing_info(subdomain);
}

void MtoNRebalancer::move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings)
{
    subdomainCreator.declare_all_subdomain_parts();
    mBulkData.modification_begin();
    change_parts_on_entities_on_all_subdomains(mappings);
    mBulkData.modification_end();
}

void MtoNRebalancer::change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping)
{
    for(int i = 0; i < subdomainCreator.get_num_subdomains(); i++)
    {
        if(subdomain_proc_mapping[i] == static_cast<unsigned>(mBulkData.parallel_rank()))
        {
            std::vector<stk::mesh::Entity> entities = get_entities_for_subdomain(i);
            subdomainCreator.move_entities_into_subdomain_part(i, entities);
        }
    }
}

stk::mesh::BulkData& MtoNRebalancer::get_bulk()
{
    return mBulkData;
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
    subdomainCreator.create_subdomain_and_write(filename,
                                                subdomain,
                                                global_num_nodes,
                                                global_num_elems,
                                                nodeSharingInfo,
                                                numSteps,
                                                timeStep);
}

}}}

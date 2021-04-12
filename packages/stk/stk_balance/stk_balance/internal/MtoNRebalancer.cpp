#include "MtoNRebalancer.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/entityDataToField.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_balance/internal/M2NDecomposer.hpp>
#include <stk_balance/setup/M2NParser.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

namespace stk {
namespace balance {
namespace internal {

MtoNRebalancer::MtoNRebalancer(stk::mesh::BulkData &bulkData,
                               stk::mesh::Field<double> &targetField,
                               M2NDecomposer &decomposer,
                               const stk::balance::M2NParsedOptions &parsedOptions)
  : m_bulkData(bulkData),
    m_subdomainCreator(bulkData, parsedOptions.targetNumProcs),
    m_decomposer(decomposer),
    m_targetDecompField(targetField),
    m_parsedOptions(parsedOptions)
{
}

MtoNRebalancer::MtoNRebalancer(stk::io::StkMeshIoBroker& ioBroker,
                               stk::mesh::Field<double> &targetField,
                               M2NDecomposer &decomposer,
                               const stk::balance::M2NParsedOptions &parsedOptions)
  : m_bulkData(ioBroker.bulk_data()),
    m_subdomainCreator(ioBroker, parsedOptions.targetNumProcs),
    m_decomposer(decomposer),
    m_targetDecompField(targetField),
    m_parsedOptions(parsedOptions)
{
}

MtoNRebalancer::~MtoNRebalancer() {}

void
MtoNRebalancer::decompose_mesh()
{
  m_decomp = m_decomposer.get_partition();
}

std::vector<unsigned>
MtoNRebalancer::map_new_subdomains_to_original_processors()
{
  return m_decomposer.map_new_subdomains_to_original_processors();
}

void MtoNRebalancer::move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc)
{
    store_off_target_proc_on_elements_before_moving_subdomains();
    stk::balance::internal::rebalance(m_bulkData, target_proc_to_starting_proc, m_decomp);
    move_entities_into_mapped_subdomain_parts(target_proc_to_starting_proc);
}

stk::io::EntitySharingInfo MtoNRebalancer::get_node_sharing_info(unsigned subdomain)
{
    return m_subdomainCreator.get_node_sharing_info(subdomain);
}

void MtoNRebalancer::move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings)
{
    m_subdomainCreator.declare_all_final_subdomain_parts();
    m_bulkData.modification_begin();
    change_parts_on_entities_on_all_subdomains(mappings);
    m_bulkData.modification_end();
}

void MtoNRebalancer::change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping)
{
    for(int i = 0; i < m_subdomainCreator.get_num_final_subdomains(); i++)
    {
        if(subdomain_proc_mapping[i] == static_cast<unsigned>(m_bulkData.parallel_rank()))
        {
            std::vector<stk::mesh::Entity> entities = get_entities_for_subdomain(i);
            m_subdomainCreator.move_entities_into_final_subdomain_part(i, entities);
        }
    }
}

stk::mesh::BulkData& MtoNRebalancer::get_bulk()
{
    return m_bulkData;
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
    double *bucketSubdomainData = static_cast<double*>(stk::mesh::field_data(m_targetDecompField, bucket));
    for(size_t k = 0; k < bucket.size(); k++)
    {
        if(bucketSubdomainData[k] == static_cast<double>(subdomain_num))
            entities.push_back(bucket[k]);
    }
}

stk::mesh::MetaData& MtoNRebalancer::get_meta()
{
    return m_bulkData.mesh_meta_data();
}

void MtoNRebalancer::store_off_target_proc_on_elements_before_moving_subdomains()
{
    stk::balance::internal::put_entity_data_to_field(m_decomp, &m_targetDecompField);
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
    m_subdomainCreator.create_subdomain_and_write(filename,
                                                subdomain,
                                                global_num_nodes,
                                                global_num_elems,
                                                nodeSharingInfo,
                                                numSteps,
                                                timeStep);
}

int MtoNRebalancer::get_num_target_subdomains()
{
  return m_subdomainCreator.get_num_final_subdomains();
}

}}}

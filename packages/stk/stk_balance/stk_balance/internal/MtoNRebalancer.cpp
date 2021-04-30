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

void MtoNRebalancer::move_final_subdomains_onto_this_processor(const std::vector<unsigned>& originalProcessorForEachFinalSubdomain)
{
    store_off_target_proc_on_elements_before_moving_subdomains();
    move_entities_into_mapped_subdomain_parts(originalProcessorForEachFinalSubdomain);
    stk::balance::internal::rebalance(m_bulkData, originalProcessorForEachFinalSubdomain, m_decomp);
}

std::vector<unsigned>
MtoNRebalancer::get_final_subdomains_for_this_processor()
{
  const std::vector<unsigned> originalProcessorsForEachNewSubdomain = map_new_subdomains_to_original_processors();

  std::vector<unsigned> finalSubdomainsForThisProc;
  for (unsigned i = 0; i < originalProcessorsForEachNewSubdomain.size(); ++i) {
    if (originalProcessorsForEachNewSubdomain[i] == static_cast<unsigned>(m_bulkData.parallel_rank())) {
      finalSubdomainsForThisProc.push_back(i);
    }
  }

  finalSubdomainsForThisProc.resize(m_decomposer.num_required_subdomains_for_each_proc(), SubdomainCreator::INVALID_SUBDOMAIN);

  return finalSubdomainsForThisProc;
}

void MtoNRebalancer::move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings)
{
  const stk::mesh::PartVector subdomainParts = m_subdomainCreator.declare_all_final_subdomain_parts();

  m_bulkData.modification_begin();
  for (const stk::mesh::EntityProc & entityProc : m_decomp) {
    m_bulkData.change_entity_parts(entityProc.first, stk::mesh::PartVector{subdomainParts[entityProc.second]});
  }
  m_bulkData.modification_end();
}

stk::mesh::BulkData& MtoNRebalancer::get_bulk()
{
    return m_bulkData;
}

SubdomainCreator&
MtoNRebalancer::get_subdomain_creator()
{
   return m_subdomainCreator;
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

void MtoNRebalancer::create_subdomain_and_write(const std::string &filename, unsigned subdomain,
                                                int global_num_nodes, int global_num_elems,
                                                int numSteps, double timeStep)
{
    m_subdomainCreator.create_subdomain_and_write(filename, subdomain, global_num_nodes, global_num_elems, numSteps, timeStep);
}

int MtoNRebalancer::get_num_target_subdomains()
{
  return m_subdomainCreator.get_num_final_subdomains();
}

}}}

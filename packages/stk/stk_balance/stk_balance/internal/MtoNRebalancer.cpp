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
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <memory>

namespace stk {
namespace balance {
namespace internal {

MtoNRebalancer::MtoNRebalancer(stk::io::StkMeshIoBroker& ioBroker,
                               stk::mesh::Field<unsigned> &targetField,
                               M2NDecomposer &decomposer,
                               const stk::balance::M2NParsedOptions &parsedOptions)
  : m_bulkData(ioBroker.bulk_data()),
    m_subdomainCreator(ioBroker, parsedOptions.targetNumProcs),
    m_decomposer(decomposer),
    m_inputMeshTargetDecompField(targetField),
    m_parsedOptions(parsedOptions)
{
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(m_bulkData, counts);
  m_globalNumNodes = counts[stk::topology::NODE_RANK];
  m_globalNumElems = counts[stk::topology::ELEM_RANK];
}

void
MtoNRebalancer::rebalance(int numSteps, double timeStep)
{
  decompose_mesh();
  map_new_subdomains_to_original_processors();
  store_final_decomp_on_elements();

  const std::vector<std::vector<unsigned>> targetSubdomainsForEachBatch = get_target_subdomains_for_each_batch();

  for (const std::vector<unsigned> & targetSubdomains : targetSubdomainsForEachBatch) {
    OutputMesh outputMesh = clone_target_subdomains(targetSubdomains);
    move_subdomain_to_owning_processor(outputMesh);
    create_subdomain_and_write(m_parsedOptions.inFile, targetSubdomains, outputMesh, numSteps, timeStep);
  }
}

void
MtoNRebalancer::decompose_mesh()
{
  m_decomp = m_decomposer.get_partition();
}

void
MtoNRebalancer::map_new_subdomains_to_original_processors()
{
  m_ownerForEachFinalSubdomain = m_decomposer.map_new_subdomains_to_original_processors();
}

void MtoNRebalancer::store_final_decomp_on_elements()
{
    store_off_target_proc_on_elements_before_moving_subdomains();
    move_entities_into_mapped_subdomain_parts();
}

std::vector<std::vector<unsigned>>
MtoNRebalancer::get_target_subdomains_for_each_batch()
{
  const unsigned numberOfBatchesToProcess = m_decomposer.num_required_subdomains_for_each_proc();

  std::vector<std::vector<unsigned>> targetSubdomainsForEachBatch(numberOfBatchesToProcess);
  for (int proc = 0; proc < m_bulkData.parallel_size(); ++proc) {
    unsigned batch = 0;
    for (unsigned i = 0; i < m_ownerForEachFinalSubdomain.size(); ++i) {
      if (m_ownerForEachFinalSubdomain[i] == static_cast<unsigned>(proc)) {
        targetSubdomainsForEachBatch[batch++].push_back(i);
      }
    }
  }

  for (unsigned batch = 0; batch < numberOfBatchesToProcess; ++batch) {
    targetSubdomainsForEachBatch[batch].resize(m_bulkData.parallel_size(), SubdomainCreator::INVALID_SUBDOMAIN);
  }

  return targetSubdomainsForEachBatch;
}

OutputMesh
MtoNRebalancer::clone_target_subdomains(const std::vector<unsigned> & targetSubdomains)
{
  OutputMesh outputMesh(m_bulkData.parallel());

  stk::mesh::Selector subdomainSelector;
  for (unsigned subdomain: targetSubdomains) {
    if (stk::transfer_utils::is_valid_subdomain(subdomain)) {
      subdomainSelector |= *m_bulkData.mesh_meta_data().get_part(m_subdomainCreator.get_subdomain_part_name(subdomain));
    }
  }

  stk::tools::copy_mesh(m_bulkData, subdomainSelector, outputMesh.bulk());

  return outputMesh;
}

void
MtoNRebalancer::move_subdomain_to_owning_processor(OutputMesh & outputMesh)
{
  stk::mesh::Field<unsigned> & outputMeshTargetDecompField =
      *reinterpret_cast<stk::mesh::Field<unsigned> *>(outputMesh.meta().get_field(stk::topology::ELEM_RANK,
                                                                                  m_inputMeshTargetDecompField.name()));

  stk::mesh::EntityProcVec subdomainDecomp;
  const stk::mesh::Part & locallyOwned = outputMesh.meta().locally_owned_part();
  for (const stk::mesh::Bucket * bucket : outputMesh.bulk().get_buckets(stk::topology::ELEM_RANK, locallyOwned)) {
    for (const stk::mesh::Entity & elem : *bucket) {
      const unsigned * destination = stk::mesh::field_data(outputMeshTargetDecompField, elem);
      subdomainDecomp.emplace_back(elem, *destination);
    }
  }

  stk::balance::internal::rebalance(outputMesh.bulk(), m_ownerForEachFinalSubdomain, subdomainDecomp);
}

void MtoNRebalancer::move_entities_into_mapped_subdomain_parts()
{
  const stk::mesh::PartVector & subdomainParts = m_subdomainCreator.declare_all_final_subdomain_parts();

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

stk::mesh::MetaData& MtoNRebalancer::get_meta()
{
    return m_bulkData.mesh_meta_data();
}

void MtoNRebalancer::store_off_target_proc_on_elements_before_moving_subdomains()
{
    stk::balance::internal::put_entity_data_to_field(m_decomp, &m_inputMeshTargetDecompField);
}

void MtoNRebalancer::create_subdomain_and_write(const std::string &filename, const std::vector<unsigned> & targetSubdomains,
                                                OutputMesh & outputMesh, int numSteps, double timeStep)
{
    m_subdomainCreator.create_subdomain_and_write(filename, targetSubdomains, m_ownerForEachFinalSubdomain,
                                                  outputMesh, m_globalNumNodes, m_globalNumElems, numSteps, timeStep);
}

int MtoNRebalancer::get_num_target_subdomains()
{
  return m_subdomainCreator.get_num_final_subdomains();
}

}}}

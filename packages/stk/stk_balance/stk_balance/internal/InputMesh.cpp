// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "InputMesh.hpp"
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_balance/internal/Decomposer.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_balance/internal/RebalanceTransientFieldTransferById.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <numeric>

namespace stk {
namespace balance {

Decomposer * make_decomposer(stk::mesh::BulkData& bulkData,
                             const stk::balance::BalanceSettings& balanceSettings)
{
  if (balanceSettings.get_use_nested_decomp()) {
    return new NestedDecomposer(bulkData, balanceSettings);
  }
  else {
    return new DefaultDecomposer(bulkData, balanceSettings);
  }
}

InputMesh::InputMesh(stk::io::StkMeshIoBroker& ioBroker,
                     const stk::balance::BalanceSettings &balanceSettings)
  : m_ioBroker(ioBroker),
    m_bulk(ioBroker.bulk_data()),
    m_meta(ioBroker.bulk_data().mesh_meta_data()),
    m_balanceSettings(balanceSettings),
    m_decomp(ioBroker.bulk_data())
{
  m_decomposer = make_decomposer(m_bulk, balanceSettings);

  initialize_mesh_counts();
  compute_output_partition();
  compute_owner_for_each_output_subdomain();
  declare_all_output_subdomain_parts();
  move_elements_into_output_subdomain_parts();
}

InputMesh::~InputMesh()
{
  delete m_decomposer;
}

unsigned
InputMesh::get_num_output_processors() const
{
  return m_balanceSettings.get_num_output_processors();
}

std::string
InputMesh::get_output_file_name() const
{
  return m_balanceSettings.get_output_filename();
}

void InputMesh::initialize_mesh_counts()
{
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(m_bulk, counts);
  m_globalNumNodes = counts[stk::topology::NODE_RANK];
  m_globalNumElems = counts[stk::topology::ELEM_RANK];
}

void
InputMesh::compute_output_partition()
{
  m_decomp = m_decomposer->get_partition();
}

void
InputMesh::compute_owner_for_each_output_subdomain()
{
  m_ownerForEachFinalSubdomain = m_decomposer->map_new_subdomains_to_original_processors();
}

bool
InputMesh::is_connected_to_element(stk::mesh::Entity entity)
{
  unsigned numElements = m_bulk.num_elements(entity);
  const stk::mesh::Entity *elems = m_bulk.begin_elements(entity);
  for (unsigned i = 0; i < numElements; ++i) {
    if (m_bulk.bucket(elems[i]).owned()) {
      return true;
    }
  }
  return false;
}

stk::mesh::EntityVector
InputMesh::get_all_orphans()
{
  stk::mesh::EntityVector orphanedEntities;

  for (stk::mesh::EntityRank rank : {stk::topology::FACE_RANK, stk::topology::EDGE_RANK, stk::topology::NODE_RANK}) {
    for (const stk::mesh::Bucket* bucket : m_bulk.buckets(rank)) {
      if (bucket->owned()) {
        for (stk::mesh::Entity entity : *bucket) {
          if (not is_connected_to_element(entity)) {
            orphanedEntities.push_back(entity);
          }
        }
      }
    }
  }

  return orphanedEntities;
}

void
InputMesh::move_elements_into_output_subdomain_parts()
{
  std::map<stk::mesh::Entity, int> finalOwner;

  const stk::mesh::BucketVector & ownedElementBuckets = m_bulk.get_buckets(stk::topology::ELEM_RANK,
                                                                           m_bulk.mesh_meta_data().locally_owned_part());
  for (const stk::mesh::Bucket * bucket : ownedElementBuckets) {
    for (const stk::mesh::Entity & elem : *bucket) {
      finalOwner[elem] = m_bulk.parallel_owner_rank(elem);
    }
  }

  const stk::mesh::EntityProcVec newDecomp = m_decomp.get_decomposition();
  for (const stk::mesh::EntityProc & entityProc : newDecomp) {
    finalOwner[entityProc.first] = entityProc.second;
  }

  const stk::mesh::EntityVector orphanedEntities = get_all_orphans();
  for (const stk::mesh::Entity & entity : orphanedEntities) {
    finalOwner[entity] = 0;
  }

  stk::mesh::EntityVector entities(finalOwner.size());
  std::vector<stk::mesh::PartVector> addParts(finalOwner.size(), stk::mesh::PartVector{nullptr});
  std::vector<stk::mesh::PartVector> removeParts(finalOwner.size());

  unsigned entityIdx = 0;
  for (const auto & entityProc : finalOwner) {
    entities[entityIdx] = entityProc.first;
    addParts[entityIdx][0] = m_subdomainParts[entityProc.second];
    ++entityIdx;
  }

  m_bulk.batch_change_entity_parts(entities, addParts, removeParts);
}

void
InputMesh::declare_all_output_subdomain_parts()
{
  const unsigned numOutputSubdomains = m_balanceSettings.get_num_output_processors();
  m_subdomainParts.clear();
  for (unsigned i = 0; i < numOutputSubdomains; ++i) {
    std::string partNameForSubdomain = get_subdomain_part_name(i);
    m_subdomainParts.push_back(&m_meta.declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK));
  }
}

std::string
InputMesh::get_subdomain_part_name(unsigned subdomainId) const
{
  return "subdomain_" + std::to_string(subdomainId);
}

std::vector<std::vector<unsigned>>
InputMesh::get_output_subdomains_for_each_batch() const
{
  const unsigned numberOfBatchesToProcess = m_decomposer->num_required_subdomains_for_each_proc();

  std::vector<std::vector<unsigned>> targetSubdomainsForEachBatch(numberOfBatchesToProcess);
  for (int proc = 0; proc < m_bulk.parallel_size(); ++proc) {
    unsigned batch = 0;
    for (unsigned i = 0; i < m_ownerForEachFinalSubdomain.size(); ++i) {
      if (m_ownerForEachFinalSubdomain[i] == static_cast<unsigned>(proc)) {
        targetSubdomainsForEachBatch[batch++].push_back(i);
      }
    }
  }

  for (unsigned batch = 0; batch < numberOfBatchesToProcess; ++batch) {
    targetSubdomainsForEachBatch[batch].resize(m_bulk.parallel_size(), INVALID_SUBDOMAIN);
  }

  return targetSubdomainsForEachBatch;
}

struct SubdomainNodeSharing {
  std::vector<stk::mesh::Entity> nodes;
  std::vector<unsigned> subdomain;
  unsigned size() const { return nodes.size(); }
  void addSharing(const stk::mesh::EntityVector & sharedNodes, unsigned otherSubdomain) {
    nodes.insert(nodes.end(), sharedNodes.begin(), sharedNodes.end());
    subdomain.resize(nodes.size(), otherSubdomain);
  }
};

stk::io::EntitySharingInfo
InputMesh::get_node_sharing_info(unsigned /*mySubdomain*/, const std::vector<unsigned>& subdomainsInBatch) const
{
  stk::CommSparse commSparse(m_bulk.parallel());

  auto isLocalSubdomain = [&](unsigned subdomain) {
    if (is_valid_subdomain(subdomain)) {
      const stk::mesh::Part* subdomainPart = m_subdomainParts[subdomain];
      const stk::mesh::BucketVector& subdomainNodeBuckets = m_bulk.get_buckets(stk::topology::NODE_RANK, *subdomainPart);
      return !subdomainNodeBuckets.empty();
    }
    return false;
  };

  std::vector<unsigned> localSubdomainsInBatch;
  std::copy_if(subdomainsInBatch.begin(), subdomainsInBatch.end(), std::back_inserter(localSubdomainsInBatch), isLocalSubdomain);

  const unsigned numOutputSubdomains = m_balanceSettings.get_num_output_processors();
  std::vector<unsigned> localSubdomains;
  std::vector<unsigned> allSubdomains(numOutputSubdomains);
  std::iota(allSubdomains.begin(), allSubdomains.end(), 0);
  std::copy_if(allSubdomains.begin(), allSubdomains.end(), std::back_inserter(localSubdomains), isLocalSubdomain);

  std::vector<SubdomainNodeSharing> nodeSharingForEachSubdomain(numOutputSubdomains);
  for (unsigned targetSubdomain : localSubdomainsInBatch) {
    for (unsigned otherSubdomain : localSubdomains) {
      if (targetSubdomain != otherSubdomain) {
        const stk::mesh::EntityVector sharedNodesThisSubdomain = get_nodes_shared_between_subdomains(targetSubdomain,
                                                                                                     otherSubdomain);
        nodeSharingForEachSubdomain[targetSubdomain].addSharing(sharedNodesThisSubdomain, otherSubdomain);
      }
    }
  }

  pack_and_communicate(commSparse,
    [&]()
    {
      for (unsigned targetSubdomain : localSubdomainsInBatch) {
        unsigned destinationProc = m_ownerForEachFinalSubdomain[targetSubdomain];
        stk::CommBuffer& buf = commSparse.send_buffer(destinationProc);
        for (unsigned i = 0; i < nodeSharingForEachSubdomain[targetSubdomain].size(); ++i) {
          buf.pack<stk::mesh::EntityId>(m_bulk.identifier(nodeSharingForEachSubdomain[targetSubdomain].nodes[i]));
          buf.pack<unsigned>(nodeSharingForEachSubdomain[targetSubdomain].subdomain[i]);
        }
      }
    });

  stk::io::EntitySharingInfo nodeSharingInfo;
  for(int proc = 0; proc < commSparse.parallel_size(); ++proc) {
    while (commSparse.recv_buffer(proc).remaining()) {
      stk::mesh::EntityId nodeId;
      unsigned sharingSubdomain = 0;
      commSparse.recv_buffer(proc).unpack<stk::mesh::EntityId>(nodeId);
      commSparse.recv_buffer(proc).unpack<unsigned>(sharingSubdomain);
      nodeSharingInfo.push_back({nodeId, sharingSubdomain});
    }
  }

  stk::util::sort_and_unique(nodeSharingInfo);

  return nodeSharingInfo;
}

stk::mesh::EntityVector
InputMesh::get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index) const
{
    int subdomain1 = std::min(this_subdomain_index, other_subdomain_index);
    int subdomain2 = std::max(this_subdomain_index, other_subdomain_index);
    stk::mesh::Selector subdomainIntersection = *m_subdomainParts[subdomain1] & *m_subdomainParts[subdomain2];
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(m_bulk, stk::topology::NODE_RANK, subdomainIntersection, nodes);
    return nodes;
}

}
}

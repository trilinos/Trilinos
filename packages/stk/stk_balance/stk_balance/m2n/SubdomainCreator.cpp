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
#include "SubdomainCreator.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <numeric>

namespace stk {
namespace balance {
namespace m2n {

constexpr unsigned SubdomainCreator::INVALID_SUBDOMAIN;

SubdomainCreator::SubdomainCreator(stk::mesh::BulkData &bulk, int numTarget)
  : m_inputMeta(bulk.mesh_meta_data()),
    m_inputBulk(bulk),
    m_numFinalSubdomains(numTarget),
    m_transientIo(nullptr)
{
}

SubdomainCreator::SubdomainCreator(stk::io::StkMeshIoBroker& ioBroker, int numTarget)
  : m_inputMeta(ioBroker.meta_data()),
    m_inputBulk(ioBroker.bulk_data()),
    m_numFinalSubdomains(numTarget),
    m_transientIo(new stk::transfer_utils::MtoNTransientFieldTransferById(ioBroker, numTarget))
{
}

SubdomainCreator::~SubdomainCreator() 
{
    delete m_transientIo;
}

std::string SubdomainCreator::get_subdomain_part_name(int subdomainId)
{
  return "subdomain_" + std::to_string(subdomainId);
}

const stk::mesh::PartVector & SubdomainCreator::declare_all_final_subdomain_parts()
{
  m_subdomainParts.clear();
  for (unsigned i = 0; i < m_numFinalSubdomains; ++i) {
    std::string partNameForSubdomain = get_subdomain_part_name(i);
    m_subdomainParts.push_back(&m_inputMeta.declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK));
  }
  return m_subdomainParts;
}

stk::mesh::Part* SubdomainCreator::get_subdomain_part(size_t subdomain)
{
  ThrowAssert(subdomain < m_subdomainParts.size());
  return m_subdomainParts[subdomain];
}

stk::mesh::PartVector SubdomainCreator::get_parts_to_add_for_subdomain(size_t subdomain_num)
{
    stk::mesh::Part* part = get_subdomain_part(subdomain_num);
    return {part};
}

stk::mesh::EntityVector SubdomainCreator::get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index)
{
    int subdomain1 = std::min(this_subdomain_index, other_subdomain_index);
    int subdomain2 = std::max(this_subdomain_index, other_subdomain_index);
    stk::mesh::Selector subdomainIntersection = *get_subdomain_part(subdomain1) & *get_subdomain_part(subdomain2);
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(m_inputBulk, stk::topology::NODE_RANK, subdomainIntersection, nodes);
    return nodes;
}

void SubdomainCreator::fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num)
{
    stk::mesh::EntityVector nodes = get_nodes_shared_between_subdomains(this_subdomain_num, other_subdomain_num);
    shared_nodes.insert(shared_nodes.end(), nodes.begin(), nodes.end());
    procs_for_shared_nodes.resize(shared_nodes.size(), other_subdomain_num);
}

void SubdomainCreator::fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes)
{
  for (unsigned other_subdomain_num=0;other_subdomain_num<m_numFinalSubdomains;++other_subdomain_num) {
    if (this_subdomain_num != other_subdomain_num) {
      fill_shared_node_proc_info(shared_nodes, procs_for_shared_nodes, this_subdomain_num, other_subdomain_num);
    }
  }
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

stk::io::EntitySharingInfo SubdomainCreator::get_node_sharing_info(unsigned mySubdomain,
                                                                   const std::vector<unsigned>& subdomainsInBatch,
                                                                   const std::vector<unsigned>& ownerForEachFinalSubdomain)
{
  stk::CommSparse commSparse(m_inputBulk.parallel());

  auto isLocalSubdomain = [&](unsigned subdomain) {
    if (stk::transfer_utils::is_valid_subdomain(subdomain)) {
      const stk::mesh::Part* subdomainPart = get_subdomain_part(subdomain);
      const stk::mesh::BucketVector& subdomainNodeBuckets = m_inputBulk.get_buckets(stk::topology::NODE_RANK, *subdomainPart);
      return !subdomainNodeBuckets.empty();
    }
    return false;
  };

  std::vector<unsigned> localSubdomainsInBatch;
  std::copy_if(subdomainsInBatch.begin(), subdomainsInBatch.end(), std::back_inserter(localSubdomainsInBatch), isLocalSubdomain);

  std::vector<unsigned> localSubdomains;
  std::vector<unsigned> allSubdomains(m_numFinalSubdomains);
  std::iota(allSubdomains.begin(), allSubdomains.end(), 0);
  std::copy_if(allSubdomains.begin(), allSubdomains.end(), std::back_inserter(localSubdomains), isLocalSubdomain);

  std::vector<SubdomainNodeSharing> nodeSharingForEachSubdomain(m_numFinalSubdomains);
  for (unsigned targetSubdomain : localSubdomainsInBatch) {
    for (unsigned otherSubdomain : localSubdomains) {
      if (targetSubdomain != otherSubdomain) {
        const stk::mesh::EntityVector sharedNodesThisSubdomain = this->get_nodes_shared_between_subdomains(targetSubdomain,
                                                                                                           otherSubdomain);
        nodeSharingForEachSubdomain[targetSubdomain].addSharing(sharedNodesThisSubdomain, otherSubdomain);
      }
    }
  }

  pack_and_communicate(commSparse,
    [&]()
    {
      for (unsigned targetSubdomain : localSubdomainsInBatch) {
        unsigned destinationProc = ownerForEachFinalSubdomain[targetSubdomain];
        stk::CommBuffer& buf = commSparse.send_buffer(destinationProc);
        for (unsigned i = 0; i < nodeSharingForEachSubdomain[targetSubdomain].size(); ++i) {
          buf.pack<stk::mesh::EntityId>(this->m_inputBulk.identifier(nodeSharingForEachSubdomain[targetSubdomain].nodes[i]));
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

void SubdomainCreator::create_subdomain_and_write(const std::string &filename, const std::vector<unsigned> & targetSubdomains,
                                                  const std::vector<unsigned> & ownerForEachFinalSubdomain,
                                                  OutputMesh & outputMesh, unsigned globalNumNodes, unsigned globalNumElems,
                                                  int numSteps, double timeStep)
{
  const unsigned mySubdomain = targetSubdomains[m_inputBulk.parallel_rank()];
  stk::io::EntitySharingInfo nodeSharingInfo = get_node_sharing_info(mySubdomain, targetSubdomains, ownerForEachFinalSubdomain);
  outputMesh.bulk().switch_to_serial_mesh();

  if (m_transientIo == nullptr) {
    if (stk::transfer_utils::is_valid_subdomain(mySubdomain)) {
      stk::io::write_file_for_subdomain(filename, mySubdomain, m_numFinalSubdomains, globalNumNodes, globalNumElems,
                                        outputMesh.bulk(), nodeSharingInfo, numSteps, timeStep);
    }
  }
  else {
    m_transientIo->setup_subdomain(outputMesh.bulk(), filename, mySubdomain, nodeSharingInfo, globalNumNodes, globalNumElems);
    m_transientIo->transfer_and_write_transient_data(mySubdomain);
  }
}     

}}}

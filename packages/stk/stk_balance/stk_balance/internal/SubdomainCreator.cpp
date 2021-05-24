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

namespace stk {
namespace balance {
namespace internal {

constexpr unsigned SubdomainCreator::INVALID_SUBDOMAIN;

SubdomainCreator::SubdomainCreator(stk::mesh::BulkData &bulk, int numTarget)
  : mMeta(bulk.mesh_meta_data()),
    mBulk(bulk),
    mNumFinalSubdomains(numTarget),
    mTransientIo(nullptr)
{
}

SubdomainCreator::SubdomainCreator(stk::io::StkMeshIoBroker& ioBroker, int numTarget)
  : mMeta(ioBroker.meta_data()),
    mBulk(ioBroker.bulk_data()),
    mNumFinalSubdomains(numTarget),
    mTransientIo(new stk::transfer_utils::MtoNTransientFieldTransferById(ioBroker, numTarget))
{
}

SubdomainCreator::~SubdomainCreator() 
{
    delete mTransientIo;
}

std::string SubdomainCreator::getSubdomainPartName(int subdomainId)
{
    std::ostringstream out;
    out << "subdomain_" << subdomainId;
    return out.str();
}

stk::mesh::PartVector SubdomainCreator::declare_all_final_subdomain_parts()
{
  stk::mesh::PartVector subdomainParts;
  for (int i = 0; i < mNumFinalSubdomains; ++i) {
    std::string partNameForSubdomain = getSubdomainPartName(i);
    subdomainParts.push_back(&mMeta.declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK));
  }
  return subdomainParts;
}

stk::mesh::Part* SubdomainCreator::get_subdomain_part(size_t subdomain_num)
{
    std::string partNameForSubdomain = getSubdomainPartName(subdomain_num);
    return mMeta.get_part(partNameForSubdomain);
}

void SubdomainCreator::move_entities_into_final_subdomain_part(size_t i, const stk::mesh::EntityVector &entities)
{
    stk::mesh::PartVector partVector = get_parts_to_add_for_subdomain(i);
    for(size_t j = 0; j < entities.size(); j++)
        mBulk.change_entity_parts(entities[j], partVector);
}

stk::mesh::PartVector SubdomainCreator::get_parts_to_add_for_subdomain(size_t subdomain_num)
{
    stk::mesh::Part* part = get_subdomain_part(subdomain_num);
    return {part};
}

stk::io::EntitySharingInfo SubdomainCreator::get_node_sharing_info(unsigned subdomain)
{
  stk::io::EntitySharingInfo nodeSharingInfo;
  if (stk::transfer_utils::is_valid_subdomain(subdomain)) {
    stk::mesh::EntityVector sharedNodes;
    std::vector<int> procsForSharedNodes;
    fill_shared_node_info_for_this_subdomain(subdomain, sharedNodes, procsForSharedNodes);

    for(size_t nodeIndex = 0; nodeIndex < sharedNodes.size(); nodeIndex++)
    {
      stk::mesh::EntityId nodeId = mBulk.identifier(sharedNodes[nodeIndex]);
      nodeSharingInfo.push_back({nodeId, procsForSharedNodes[nodeIndex]});
    }
  }
  return nodeSharingInfo;
}

stk::mesh::EntityVector SubdomainCreator::get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index)
{
    stk::mesh::Selector selected_nodes = *get_subdomain_part(this_subdomain_index) &
                                         *get_subdomain_part(other_subdomain_index);
    const bool sortById = true;
    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(mBulk, stk::topology::NODE_RANK, selected_nodes, nodes, sortById);
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
    for(int other_subdomain_num=0;other_subdomain_num<mNumFinalSubdomains;++other_subdomain_num)
    {
        if(static_cast<int>(this_subdomain_num) != other_subdomain_num)
            fill_shared_node_proc_info(shared_nodes, procs_for_shared_nodes, this_subdomain_num, other_subdomain_num);
    }
}

void SubdomainCreator::create_subdomain_and_write(const std::string &filename, unsigned subdomain,
                                                  int global_num_nodes, int global_num_elems,
                                                  int numSteps, double timeStep)
{
  stk::mesh::MetaData newMeta;
  stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);

  const stk::mesh::Selector subdomainSelector = stk::transfer_utils::is_valid_subdomain(subdomain) ? *mMeta.get_part(getSubdomainPartName(subdomain))
                                                                                                   : stk::mesh::Selector();
  stk::tools::copy_mesh(mBulk, subdomainSelector, newBulkData);

  stk::io::EntitySharingInfo nodeSharingInfo = get_node_sharing_info(subdomain);

  if (mTransientIo == nullptr) {
    if (stk::transfer_utils::is_valid_subdomain(subdomain)) {
      stk::io::write_file_for_subdomain(filename, subdomain, mNumFinalSubdomains, global_num_nodes, global_num_elems,
                                        newBulkData, nodeSharingInfo, numSteps, timeStep);
    }
  }
  else {
    mTransientIo->setup_subdomain(newBulkData, filename, subdomain, nodeSharingInfo, global_num_nodes, global_num_elems);
    mTransientIo->transfer_and_write_transient_data(subdomain);
  }
}     

}}}

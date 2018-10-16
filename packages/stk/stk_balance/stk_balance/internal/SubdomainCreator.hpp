// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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
#ifndef STK_BALANCE_SUBDOMAINCREATOR_HPP
#define STK_BALANCE_SUBDOMAINCREATOR_HPP

#include <string>
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <stk_io/IossBridge.hpp>

namespace stk { namespace mesh { class MetaData; }}
namespace stk { namespace mesh { class BulkData; }}

namespace stk {
namespace balance {
namespace internal {

class SubdomainCreator
{
public:
    SubdomainCreator(stk::mesh::BulkData &bulk, int numTarget) : mMeta(bulk.mesh_meta_data()), mBulk(bulk), mNumSubdomains(numTarget) {}

    void declare_all_subdomain_parts();
    void move_entities_into_subdomain_part(size_t i, const stk::mesh::EntityVector &entities);
    stk::io::EntitySharingInfo get_node_sharing_info(unsigned subdomain);
    void create_subdomain_and_write(const std::string &filename, unsigned subdomain, int global_num_nodes, int global_num_elems, const stk::io::EntitySharingInfo &nodeSharingInfo, int numSteps = -1, double timeStep = 0.0);

    int get_num_subdomains() const { return mNumSubdomains; }
private:
    std::string getSubdomainPartName(int subdomainId);
    stk::mesh::Part* get_subdomain_part(size_t subdomain_num);
    stk::mesh::PartVector get_parts_to_add_for_subdomain(size_t subdomain_num);
    stk::mesh::EntityVector get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index);
    void fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num);
    void fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes);

    stk::mesh::MetaData &mMeta;
    stk::mesh::BulkData &mBulk;
    int mNumSubdomains;
};

}}}
#endif

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
#ifndef STK_BALANCE_MTONREBALANCER_HPP
#define STK_BALANCE_MTONREBALANCER_HPP

#include "SubdomainCreator.hpp"
#include <stk_balance/balanceUtils.hpp>
#include <string>
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <stk_io/IossBridge.hpp>

namespace stk { namespace mesh { class MetaData; }}
namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace mesh { class FieldBase; }}
namespace stk { namespace mesh { class Bucket; }}

namespace stk {
namespace balance {
namespace internal {

class MtoNRebalancer
{
public:
    MtoNRebalancer(stk::mesh::BulkData &bulkData, stk::mesh::Field<double> &targetField, const stk::balance::BalanceSettings &graphSettings, int num_target_procs);
    virtual ~MtoNRebalancer();

    void generate_n_proc_decomp();
    void move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc);
    stk::io::EntitySharingInfo get_node_sharing_info(unsigned subdomain);
    void create_subdomain_and_write(const std::string &filename, unsigned subdomain, int global_num_nodes, int global_num_elems, const stk::io::EntitySharingInfo &nodeSharingInfo, int numSteps = -1, double timeStep = 0.0);
    bool does_this_proc_own_subdomain(unsigned subdomainOwner);

private:
    void move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings);
    void change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping);
    std::vector<stk::mesh::Entity> get_entities_for_subdomain(size_t subdomain_num);
    stk::mesh::EntityVector get_entitites_for_subdomain_using_field_from_buckets(size_t subdomain_num, const stk::mesh::BucketVector& buckets);
    void add_owned_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities);
    void add_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities);
    void store_off_target_proc_on_elements_before_moving_subdomains();
    stk::mesh::MetaData& get_meta();
    stk::mesh::BulkData& get_bulk();

private:
    MtoNRebalancer( const MtoNRebalancer& other );
    MtoNRebalancer& operator=( const MtoNRebalancer& other );

    stk::mesh::BulkData &mBulkData;
    SubdomainCreator subdomainCreator;
    stk::mesh::Field<double> &targetDecompField;
    stk::mesh::EntityProcVec decomp;
};

}}}
#endif

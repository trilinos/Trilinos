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
namespace stk { namespace balance { class M2NParsedOptions; }}
namespace stk { namespace balance { namespace internal { class M2NDecomposer; }}}
namespace stk { namespace io { class StkMeshIoBroker; }}

namespace stk {
namespace balance {
namespace internal {

class MtoNRebalancer
{
public:
    MtoNRebalancer(stk::io::StkMeshIoBroker& ioBroker,
                   stk::mesh::Field<unsigned> &targetField,
                   M2NDecomposer &decomposer,
                   const M2NBalanceSettings &balanceSettings);

    ~MtoNRebalancer() = default;

    void rebalance(int numSteps, double timeStep);

    void decompose_mesh();
    void map_new_subdomains_to_original_processors();
    std::vector<std::vector<unsigned>> get_target_subdomains_for_each_batch();

    void store_final_decomp_on_elements();

    OutputMesh clone_target_subdomains(const std::vector<unsigned> & targetSubdomains);

    void move_subdomain_to_owning_processor(OutputMesh & outputMesh);

    void create_subdomain_and_write(const std::string &filename, const std::vector<unsigned> & targetSubdomains,
                                    OutputMesh & outputMesh,  int numSteps = -1, double timeStep = 0.0);

    stk::mesh::MetaData& get_meta();
    stk::mesh::BulkData& get_bulk();
    SubdomainCreator& get_subdomain_creator();

    int get_num_target_subdomains();

    const std::vector<unsigned> & get_owner_for_each_final_subdomain() const { return m_ownerForEachFinalSubdomain; }

private:
    void move_entities_into_mapped_subdomain_parts();
    void store_off_target_proc_on_elements_before_moving_subdomains();

    MtoNRebalancer( const MtoNRebalancer& other );
    MtoNRebalancer& operator=( const MtoNRebalancer& other );

    stk::mesh::BulkData &m_bulkData;
    SubdomainCreator m_subdomainCreator;
    M2NDecomposer &m_decomposer;
    stk::mesh::Field<unsigned> &m_inputMeshTargetDecompField;
    const M2NBalanceSettings &m_balanceSettings;
    stk::mesh::EntityProcVec m_decomp;
    std::vector<unsigned> m_ownerForEachFinalSubdomain;
    int m_globalNumNodes;
    int m_globalNumElems;
};

}}}
#endif

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
#ifndef TRANSIENT_FIELD_TRANSFER_BY_ID_
#define TRANSIENT_FIELD_TRANSFER_BY_ID_

#include <vector>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include "stk_transfer/copy_by_id/TransferCopyById.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdStkMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/SearchByIdGeometric.hpp"

#include <stk_io/StkMeshIoBroker.hpp>

namespace Ioss { class Region; }

namespace stk {
namespace transfer_utils {

class TransientTransferByIdForRank
{
public:
    TransientTransferByIdForRank(stk::mesh::MetaData &metaA, stk::mesh::MetaData &metaB, stk::mesh::EntityRank rank);
    ~TransientTransferByIdForRank();

    void initialize();

    void do_transfer();

    stk::mesh::EntityRank get_rank() const { return mRank; }

    stk::mesh::MetaData  &get_metaA() { return mMetaA; }
    stk::mesh::MetaData  &get_metaB() { return mMetaB; }

protected:
    stk::mesh::MetaData   &mMetaA;
    stk::mesh::MetaData   &mMetaB;
    stk::mesh::EntityRank  mRank;

    stk::transfer::TransferCopyByIdStkMeshAdapter *mTransferMeshA = nullptr;
    stk::transfer::TransferCopyByIdStkMeshAdapter *mTransferMeshB = nullptr;

    stk::transfer::SearchByIdGeometric  mSearch;
    stk::transfer::TransferCopyById    *mTransfer = nullptr;

private:
    TransientTransferByIdForRank();

private:
    stk::transfer::TransferCopyByIdStkMeshAdapter *create_transfer_mesh(stk::mesh::MetaData &meta);
};

class TransientFieldTransferById
{
public:
    TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB);

    TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB, const std::vector<stk::mesh::EntityRank> &entityRanks);

    ~TransientFieldTransferById();

    void writeFields(size_t aOutFileIndex, std::vector<const stk::mesh::FieldBase *> & aTransientFields, std::vector<std::string> & aGlobalVariableNames);
    void get_field_names(size_t aOutputFileIndex, std::vector<const stk::mesh::FieldBase *> & aTransientFields, std::vector<std::string> & aGlobalVariableNames);

    size_t transfer_and_write_transient_fields(const std::string &parallelOutputMeshName,  stk::mesh::Selector & aselector);
    size_t transfer_and_write_transient_fields(const std::string &parallelOutputMeshName);

    stk::io::StkMeshIoBroker &get_brokerA() { return mBrokerA; }
    stk::io::StkMeshIoBroker &get_brokerB() { return mBrokerB; }

protected:
    stk::io::StkMeshIoBroker &mBrokerA;
    stk::io::StkMeshIoBroker &mBrokerB;
    std::vector<TransientTransferByIdForRank*> mTransfers;

private:
    TransientFieldTransferById();

    void do_transfer();

    size_t setup_output_transient_fields(const std::string &parallelOutputMeshName);

    void initialize(const std::vector<stk::mesh::EntityRank>& entityRanks);
};

struct SubDomainInfo {
    stk::mesh::BulkData* bulk = nullptr;
    std::string filename;
    stk::io::EntitySharingInfo nodeSharingInfo;
    int globalNumNodes = 0;
    int globalNumElems = 0;
    Ioss::Region* outRegion = nullptr;
    bool meshWritten = false;
    std::vector<TransientTransferByIdForRank*> transferVec;
};

class MtoNTransientFieldTransferById
{
public:
    MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain);

    MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain, 
                                   const std::vector<stk::mesh::EntityRank> &entityRanks);

    ~MtoNTransientFieldTransferById();

    void setup_subdomain(stk::mesh::BulkData& bulk, const std::string &filename, 
                         unsigned subdomain, const stk::io::EntitySharingInfo& nodeSharingInfo,
                         int global_num_nodes, int global_num_elems);

    void write_mesh_data(unsigned subdomain);
    void write_transient_data(unsigned subdomain, double timeStep);
    void transfer_transient_data(unsigned subdomain);
    void transfer_and_write_transient_data(unsigned subdomain);

private:
    MtoNTransientFieldTransferById();

    void validate_subdomain(unsigned subdomain);
    void initialize_transfer(unsigned subdomain);
    void add_qa_records(unsigned subdomain);
    void add_info_records(unsigned subdomain);
    void add_global_variables(unsigned subdomain);
    void write_global_variables(unsigned subdomain, int step);

    stk::io::StkMeshIoBroker &m_inputBroker;
    unsigned m_numSubDomain;
    std::vector<SubDomainInfo> m_subDomainInfoVec;
    std::vector<stk::mesh::EntityRank> m_entityRanks;
};

}
}

#endif // TRANSIENT_FIELD_TRANSFER_BY_ID_

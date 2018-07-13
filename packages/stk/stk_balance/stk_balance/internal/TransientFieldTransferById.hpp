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

namespace stk {
namespace balance {
namespace internal {

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

}
}
}

#endif // TRANSIENT_FIELD_TRANSFER_BY_ID_

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

#include <stk_util/registry/ProductRegistry.hpp>
#include "TransientFieldTransferById.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "internal/privateDeclarations.hpp"

namespace stk {
namespace balance {
namespace internal {

TransientTransferByIdForRank::TransientTransferByIdForRank(stk::mesh::MetaData &metaA,
                                                           stk::mesh::MetaData &metaB,
                                                           stk::mesh::EntityRank rank) :
        mMetaA(metaA),
        mMetaB(metaB),
        mRank(rank)
{
    mTransferMeshA = create_transfer_mesh(metaA);
    mTransferMeshB = create_transfer_mesh(metaB);

    mTransfer = new stk::transfer::TransferCopyById(mSearch, *mTransferMeshA, *mTransferMeshB);
}

TransientTransferByIdForRank::~TransientTransferByIdForRank()
{
    delete mTransfer;
    mTransfer = nullptr;
    delete mTransferMeshA;
    mTransferMeshA = nullptr;
    delete mTransferMeshB;
    mTransferMeshB = nullptr;
}

void TransientTransferByIdForRank::initialize()
{
    mTransfer->initialize();
}

void TransientTransferByIdForRank::do_transfer()
{
    mTransfer->apply();
}

stk::transfer::TransferCopyByIdStkMeshAdapter *TransientTransferByIdForRank::create_transfer_mesh(stk::mesh::MetaData &meta)
{
    stk::mesh::EntityVector entities;
    stk::mesh::get_selected_entities(meta.locally_owned_part(), meta.mesh_bulk_data().buckets(mRank), entities);

    stk::mesh::FieldVector fields = stk::io::get_transient_fields(meta, mRank);
    return new stk::transfer::TransferCopyByIdStkMeshAdapter(meta.mesh_bulk_data(), entities, fields);
}

TransientFieldTransferById::TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB) :
        mBrokerA(brokerA),
        mBrokerB(brokerB)
{
    std::vector<stk::mesh::EntityRank> entityRanks;

    for(const std::string &name : brokerA.meta_data().entity_rank_names())
    {
        stk::mesh::EntityRank rank = brokerA.meta_data().entity_rank(name);
        entityRanks.push_back(rank);
    }

    initialize(entityRanks);

    int dbIntSize = mBrokerA.check_integer_size_requirements();
    if(dbIntSize > 4) {
        mBrokerB.property_add(Ioss::Property("INTEGER_SIZE_API" , dbIntSize));
        mBrokerB.property_add(Ioss::Property("INTEGER_SIZE_DB" , dbIntSize));
    }
}

TransientFieldTransferById::TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA,
                                                       stk::io::StkMeshIoBroker &brokerB,
                                                       const std::vector<stk::mesh::EntityRank> &entityRanks) :
        mBrokerA(brokerA),
        mBrokerB(brokerB)
{
    initialize(entityRanks);
}

TransientFieldTransferById::~TransientFieldTransferById()
{
    for(unsigned i = 0; i < mTransfers.size(); i++)
        delete mTransfers[i];

    mTransfers.clear();
}

size_t TransientFieldTransferById::transfer_and_write_transient_fields(const std::string &parallelOutputMeshName)
{
    internal::logMessage(mBrokerA.bulk_data().parallel(), "Writing output mesh");
    size_t outputFileIndex = setup_output_transient_fields(parallelOutputMeshName);

    std::vector<stk::io::QaRecord> qaRecords = mBrokerA.get_qa_records();
    mBrokerB.add_qa_records(outputFileIndex, qaRecords);
    mBrokerB.set_name_and_version_for_qa_record(outputFileIndex, "stk_balance", stk::ProductRegistry::version());

    mBrokerB.add_info_records(outputFileIndex, mBrokerA.get_info_records());

    stk::mesh::FieldVector fieldVector = stk::io::get_transient_fields(mBrokerB.meta_data());

    std::vector<const stk::mesh::FieldBase *> transientFields(fieldVector.begin(), fieldVector.end());

    std::vector<std::string> globalVariableNames;
    mBrokerA.get_global_variable_names(globalVariableNames);
    std::vector<double> globalVariable;
    for(const std::string& globalVariableName : globalVariableNames) {
        size_t length = mBrokerA.get_global_variable_length(globalVariableName);
        mBrokerB.add_global(outputFileIndex, globalVariableName, length, Ioss::Field::DOUBLE);
    }

    if(fieldVector.empty())
    {
        mBrokerB.write_output_mesh(outputFileIndex);
    }

    for(int iStep = 0; iStep < mBrokerA.get_num_time_steps(); iStep++)
    {
        internal::logMessage(mBrokerA.bulk_data().parallel(), "Appending transient data for time step " + std::to_string(iStep));
        double readTime = mBrokerA.read_defined_input_fields_at_step(iStep + 1, nullptr);

        do_transfer();

        stk::mesh::copy_owned_to_shared(mBrokerB.bulk_data(), transientFields);

        mBrokerB.begin_output_step(outputFileIndex, readTime);
        mBrokerB.write_defined_output_fields(outputFileIndex);

        for(const std::string& globalVariableName : globalVariableNames) {
            mBrokerA.get_global(globalVariableName, globalVariable);
            mBrokerB.write_global(outputFileIndex, globalVariableName, globalVariable);
        }

        mBrokerB.end_output_step(outputFileIndex);
    }

    return outputFileIndex;
}

void TransientFieldTransferById::do_transfer()
{
    for(TransientTransferByIdForRank *transfer : mTransfers)
        transfer->do_transfer();
}

size_t TransientFieldTransferById::setup_output_transient_fields(const std::string &parallelOutputMeshName)
{
    size_t outputFileIndex = mBrokerB.create_output_mesh(parallelOutputMeshName, stk::io::WRITE_RESULTS);
    mBrokerB.set_reference_input_region(outputFileIndex, mBrokerA);

    stk::mesh::FieldVector allTransientFieldsB = stk::io::get_transient_fields(mBrokerB.meta_data());

    for(stk::mesh::FieldBase* field : allTransientFieldsB)
        mBrokerB.add_field(outputFileIndex, *field);

    return outputFileIndex;
}

void TransientFieldTransferById::initialize(const std::vector<stk::mesh::EntityRank>& entityRanks)
{
    for(stk::mesh::EntityRank rank : entityRanks)
    {
        if(stk::io::get_transient_fields(mBrokerA.meta_data(), rank).size() > 0)
        {
            TransientTransferByIdForRank *transfer = new TransientTransferByIdForRank(mBrokerA.meta_data(),
                                                                                      mBrokerB.meta_data(),
                                                                                      rank);
            transfer->initialize();
            mTransfers.push_back(transfer);
        }
    }
}

}
}
}

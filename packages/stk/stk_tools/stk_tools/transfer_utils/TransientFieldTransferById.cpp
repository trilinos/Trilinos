#include <stk_util/registry/ProductRegistry.hpp>
#include "TransientFieldTransferById.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "Ioss_Region.h"
#include "stk_io/IOHelpers.hpp"
#include <Teuchos_RCP.hpp>                         // for RCP::RCP<T>, etc
#include "Teuchos_RCPDecl.hpp"                     // for RCP

namespace stk {
namespace transfer_utils {

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

void TransientFieldTransferById::writeFields(size_t aOutFileIndex, std::vector<const stk::mesh::FieldBase *> & aTransientFields, std::vector<std::string> & aGlobalVariableNames)
{
	if(aTransientFields.empty())
	{
	    mBrokerB.write_output_mesh(aOutFileIndex);
	}

	std::vector<double> globalVariable;
	for(int iStep = 0; iStep < mBrokerA.get_num_time_steps(); iStep++)
	{
	    double readTime = mBrokerA.read_defined_input_fields_at_step(iStep + 1, nullptr);

	    do_transfer();
	    stk::mesh::copy_owned_to_shared(mBrokerB.bulk_data(), aTransientFields);

	    mBrokerB.begin_output_step(aOutFileIndex, readTime);
	    mBrokerB.write_defined_output_fields(aOutFileIndex);

	    for(const std::string& globalVariableName : aGlobalVariableNames) {
	        mBrokerA.get_global(globalVariableName, globalVariable);
	        mBrokerB.write_global(aOutFileIndex, globalVariableName, globalVariable);
	    }

	    mBrokerB.end_output_step(aOutFileIndex);
	}
}

void TransientFieldTransferById::get_field_names(size_t aOutputFileIndex, std::vector<const stk::mesh::FieldBase *> & aTransientFields, std::vector<std::string> & aGlobalVariableNames)
{
    std::vector<stk::io::QaRecord> qaRecords = mBrokerA.get_qa_records();
    mBrokerB.add_qa_records(aOutputFileIndex, qaRecords);
    mBrokerB.set_name_and_version_for_qa_record(aOutputFileIndex, "stk_balance", stk::ProductRegistry::version());

    mBrokerB.add_info_records(aOutputFileIndex, mBrokerA.get_info_records());

    stk::mesh::FieldVector fieldVector = stk::io::get_transient_fields(mBrokerB.meta_data());

    aTransientFields.assign(fieldVector.begin(), fieldVector.end());

    mBrokerA.get_global_variable_names(aGlobalVariableNames);

    for(const std::string& globalVariableName : aGlobalVariableNames) {
        size_t length = mBrokerA.get_global_variable_length(globalVariableName);
        mBrokerB.add_global(aOutputFileIndex, globalVariableName, length, Ioss::Field::DOUBLE);
    }
}

size_t TransientFieldTransferById::transfer_and_write_transient_fields(const std::string &parallelOutputMeshName)
{
    size_t outputFileIndex = setup_output_transient_fields(parallelOutputMeshName);

    std::vector<const stk::mesh::FieldBase *> transientFields;
    std::vector<std::string> globalVariableNames;

    get_field_names(outputFileIndex, transientFields, globalVariableNames);

    writeFields(outputFileIndex, transientFields, globalVariableNames);

    return outputFileIndex;
}

size_t TransientFieldTransferById::transfer_and_write_transient_fields(const std::string &parallelOutputMeshName, stk::mesh::Selector & aselector)
{
    size_t outputFileIndex = setup_output_transient_fields(parallelOutputMeshName);

    std::vector<const stk::mesh::FieldBase *> transientFields;
    std::vector<std::string> globalVariableNames;

    get_field_names(outputFileIndex, transientFields, globalVariableNames);

    mBrokerB.set_subset_selector(outputFileIndex,aselector);

    writeFields(outputFileIndex, transientFields, globalVariableNames);

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

MtoNTransientFieldTransferById::MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker& inputBroker, unsigned numSubDomain) :
        m_inputBroker(inputBroker), m_numSubDomain(numSubDomain), m_subDomainInfoVec(numSubDomain)
{
    for(const std::string &name : inputBroker.meta_data().entity_rank_names()) {
        stk::mesh::EntityRank rank = inputBroker.meta_data().entity_rank(name);
        m_entityRanks.push_back(rank);
    }
}

MtoNTransientFieldTransferById::MtoNTransientFieldTransferById(stk::io::StkMeshIoBroker &inputBroker, unsigned numSubDomain,
                                                               const std::vector<stk::mesh::EntityRank> &entityRanks) :
        m_inputBroker(inputBroker), m_numSubDomain(numSubDomain), m_subDomainInfoVec(numSubDomain), m_entityRanks(entityRanks)
{
}

MtoNTransientFieldTransferById::~MtoNTransientFieldTransferById()
{
    for(SubDomainInfo& info : m_subDomainInfoVec) {
        if(info.outRegion != nullptr) {
            stk::io::delete_selector_property(*info.outRegion);
            delete info.outRegion;
            info.outRegion = nullptr;

            for(unsigned i = 0; i < info.transferVec.size(); i++) {
                delete info.transferVec[i];
            }
            info.transferVec.clear();
        }
    }
}

void MtoNTransientFieldTransferById::initialize_transfer(unsigned subdomain)
{
    validate_subdomain(subdomain);
    for(stk::mesh::EntityRank rank : m_entityRanks) {
        if(stk::io::get_transient_fields(m_inputBroker.meta_data(), rank).size() > 0) {
            stk::mesh::MetaData& outputMeta = m_subDomainInfoVec[subdomain].bulk->mesh_meta_data();
            TransientTransferByIdForRank *transfer = new TransientTransferByIdForRank(m_inputBroker.meta_data(),
                                                                                      outputMeta,
                                                                                      rank);
            transfer->initialize();
            m_subDomainInfoVec[subdomain].transferVec.push_back(transfer);
        }
    }
}

void MtoNTransientFieldTransferById::setup_subdomain(stk::mesh::BulkData& bulk, const std::string &filename, 
                                                     unsigned subdomain, const stk::io::EntitySharingInfo& nodeSharingInfo,
                                                     int global_num_nodes, int global_num_elems)
{
    ThrowRequireMsg(subdomain < m_numSubDomain, "Invalid subdomain index: " << subdomain);
    m_subDomainInfoVec[subdomain].bulk = &bulk;
    m_subDomainInfoVec[subdomain].filename = filename;
    m_subDomainInfoVec[subdomain].nodeSharingInfo = nodeSharingInfo;
    m_subDomainInfoVec[subdomain].globalNumNodes = global_num_nodes;
    m_subDomainInfoVec[subdomain].globalNumElems = global_num_elems; 

    Ioss::DatabaseIO *dbo = stk::io::create_database_for_subdomain(filename, subdomain, m_numSubDomain);
    m_subDomainInfoVec[subdomain].outRegion = new Ioss::Region(dbo, filename);

    stk::io::add_properties_for_subdomain(bulk, *m_subDomainInfoVec[subdomain].outRegion, 
                                          subdomain, m_numSubDomain, global_num_nodes, global_num_elems);

    int dbIntSize = m_inputBroker.check_integer_size_requirements_serial();
    if(dbIntSize > 4) {
        m_subDomainInfoVec[subdomain].outRegion->property_add(Ioss::Property("INTEGER_SIZE_API", dbIntSize));
        m_subDomainInfoVec[subdomain].outRegion->property_add(Ioss::Property("INTEGER_SIZE_DB", dbIntSize));
    }
    initialize_transfer(subdomain);
}

void MtoNTransientFieldTransferById::validate_subdomain(unsigned subdomain)
{
    ThrowRequireMsg(subdomain < m_numSubDomain, "Invalid subdomain index: " << subdomain);
    ThrowRequireMsg(m_subDomainInfoVec[subdomain].bulk != nullptr, "Subdomain: " << subdomain << " has not been initialized");
}

void MtoNTransientFieldTransferById::write_mesh_data(unsigned subdomain)
{
    validate_subdomain(subdomain);
    if(!m_subDomainInfoVec[subdomain].meshWritten) {
        add_qa_records(subdomain);
        add_info_records(subdomain);
        stk::io::write_file_for_subdomain(*m_subDomainInfoVec[subdomain].outRegion, 
                                          *m_subDomainInfoVec[subdomain].bulk,
                                          m_subDomainInfoVec[subdomain].nodeSharingInfo);
        m_subDomainInfoVec[subdomain].meshWritten = true;
    }
}

void MtoNTransientFieldTransferById::write_global_variables(unsigned subdomain, int step)
{
    Ioss::Region* region = m_subDomainInfoVec[subdomain].outRegion;
    region->begin_mode(Ioss::STATE_TRANSIENT);
    region->begin_state(step);

    std::vector<std::string> globalVariableNames;
    std::vector<double> globalVariable;
    m_inputBroker.get_global_variable_names(globalVariableNames);

    Teuchos::RCP<Ioss::Region> rcpRegion;
    rcpRegion.reset(region, false);

	for(const std::string& globalVariableName : globalVariableNames) {
	    m_inputBroker.get_global(globalVariableName, globalVariable);
        stk::io::internal_write_global(rcpRegion, globalVariableName, globalVariable);
	}

    region->end_state(step);
    region->end_mode(Ioss::STATE_TRANSIENT);
}

void MtoNTransientFieldTransferById::write_transient_data(unsigned subdomain, double timeStep)
{
    validate_subdomain(subdomain);
    if(!m_subDomainInfoVec[subdomain].meshWritten) {
        write_mesh_data(subdomain);
    }

    Ioss::Region* region = m_subDomainInfoVec[subdomain].outRegion;
    int step = stk::io::write_transient_data_for_subdomain(*region,
                                                           *m_subDomainInfoVec[subdomain].bulk,
                                                           timeStep);
    write_global_variables(subdomain, step); 
}

void MtoNTransientFieldTransferById::add_qa_records(unsigned subdomain)
{
    std::vector<stk::io::QaRecord> qaRecords = m_inputBroker.get_qa_records();
    Ioss::Region* region = m_subDomainInfoVec[subdomain].outRegion;

    for(const stk::io::QaRecord &qaRec : qaRecords) {
        region->add_qa_record(qaRec.name, qaRec.version, qaRec.date, qaRec.time);
    }

    std::string codeName = "stk_balanceMtoN";
    std::string codeVersion = stk::ProductRegistry::version();
    region->property_add(Ioss::Property(std::string("code_name"), codeName));
    region->property_add(Ioss::Property(std::string("code_version"), codeVersion));
}

void MtoNTransientFieldTransferById::add_info_records(unsigned subdomain)
{
    Ioss::Region* region = m_subDomainInfoVec[subdomain].outRegion;
    region->add_information_records(m_inputBroker.get_info_records());
}

void MtoNTransientFieldTransferById::add_global_variables(unsigned subdomain)
{
    std::vector<std::string> globalVariableNames;
    m_inputBroker.get_global_variable_names(globalVariableNames);
    Ioss::Region* region = m_subDomainInfoVec[subdomain].outRegion;
    Teuchos::RCP<Ioss::Region> rcpRegion;
    rcpRegion.reset(region, false);

    for(const std::string& globalVariableName : globalVariableNames) {
        size_t length = m_inputBroker.get_global_variable_length(globalVariableName);
        stk::io::internal_add_global(rcpRegion, globalVariableName, length, Ioss::Field::DOUBLE);
    }
}

void MtoNTransientFieldTransferById::transfer_transient_data(unsigned subdomain)
{
    validate_subdomain(subdomain);
    stk::mesh::BulkData& bulk = *m_subDomainInfoVec[subdomain].bulk;
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();

    for(TransientTransferByIdForRank *transfer : m_subDomainInfoVec[subdomain].transferVec) {
       transfer->do_transfer();
    }

    stk::mesh::FieldVector fieldVector = stk::io::get_transient_fields(meta);
    std::vector<const stk::mesh::FieldBase*> transientFields;
    transientFields.assign(fieldVector.begin(), fieldVector.end());

	stk::mesh::copy_owned_to_shared(bulk, transientFields);
}

void MtoNTransientFieldTransferById::transfer_and_write_transient_data(unsigned subdomain)
{
    validate_subdomain(subdomain);
    write_mesh_data(subdomain);
    add_global_variables(subdomain);

    std::vector<double> inputTimeSteps = m_inputBroker.get_time_steps();

    for(double time : inputTimeSteps) {
        m_inputBroker.read_defined_input_fields(time);
        transfer_transient_data(subdomain);
        write_transient_data(subdomain, time);
    } 
}

}
}
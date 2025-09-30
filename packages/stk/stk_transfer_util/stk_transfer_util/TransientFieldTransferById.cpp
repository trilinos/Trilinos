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

#include <stk_util/registry/ProductRegistry.hpp>
#include "TransientFieldTransferById.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "Ioss_Region.h"
#include "stk_io/IOHelpers.hpp"

namespace stk {
namespace transfer_utils {

TransientTransferByIdForRank::TransientTransferByIdForRank(stk::mesh::MetaData &metaA,
                                                           stk::mesh::MetaData &metaB,
                                                           stk::mesh::EntityRank rank)
  : mMetaA(metaA),
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

RepeatedTransferCopyByIdStkMeshAdapter *TransientTransferByIdForRank::create_transfer_mesh(stk::mesh::MetaData &meta)
{
  stk::mesh::EntityVector entities;
  const bool sortById = true;
  stk::mesh::get_entities(meta.mesh_bulk_data(), mRank, meta.locally_owned_part(), entities, sortById);

  stk::mesh::FieldVector fields = stk::io::get_transient_fields(meta, mRank);
  return new RepeatedTransferCopyByIdStkMeshAdapter(meta.mesh_bulk_data(), entities, fields);
}

TransientFieldTransferById::TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB)
  : mBrokerA(brokerA),
    mBrokerB(brokerB)
{
  std::vector<stk::mesh::EntityRank> entityRanks;

  for (const std::string &name : brokerA.meta_data().entity_rank_names()) {
    stk::mesh::EntityRank rank = brokerA.meta_data().entity_rank(name);
    entityRanks.push_back(rank);
  }

  initialize(entityRanks);

  int dbIntSize = mBrokerA.check_integer_size_requirements();
  if (dbIntSize > 4) {
    mBrokerB.property_add(Ioss::Property("INTEGER_SIZE_API" , dbIntSize));
    mBrokerB.property_add(Ioss::Property("INTEGER_SIZE_DB" , dbIntSize));
  }
}

TransientFieldTransferById::~TransientFieldTransferById()
{
  for (unsigned i = 0; i < mTransfers.size(); i++) {
    delete mTransfers[i];
  }

  mTransfers.clear();
}

void TransientFieldTransferById::writeFields(size_t aOutFileIndex,
                                             std::vector<const stk::mesh::FieldBase *> & aTransientFields,
                                             std::vector<std::string> & aGlobalVariableNames)
{
  if (aTransientFields.empty()) {
    mBrokerB.write_output_mesh(aOutFileIndex);
  }

  std::vector<double> globalVariable;
  for (int iStep = 0; iStep < mBrokerA.get_num_time_steps(); iStep++) {
    double readTime = mBrokerA.read_defined_input_fields_at_step(iStep + 1, nullptr);

    do_transfer();
    stk::mesh::copy_owned_to_shared(mBrokerB.bulk_data(), aTransientFields);

    mBrokerB.begin_output_step(aOutFileIndex, readTime);
    mBrokerB.write_defined_output_fields(aOutFileIndex);

    for (const std::string& globalVariableName : aGlobalVariableNames) {
      mBrokerA.get_global(globalVariableName, globalVariable);
      mBrokerB.write_global(aOutFileIndex, globalVariableName, globalVariable);
    }

    mBrokerB.end_output_step(aOutFileIndex);
  }
}

void TransientFieldTransferById::get_field_names(size_t aOutputFileIndex,
                                                 std::vector<const stk::mesh::FieldBase *> & aTransientFields,
                                                 std::vector<std::string> & aGlobalVariableNames)
{
  std::vector<stk::io::QaRecord> qaRecords = mBrokerA.get_qa_records();
  mBrokerB.add_qa_records(aOutputFileIndex, qaRecords);
  mBrokerB.set_name_and_version_for_qa_record(aOutputFileIndex, "stk_balance", stk::ProductRegistry::version());

  mBrokerB.add_info_records(aOutputFileIndex, mBrokerA.get_info_records());

  stk::mesh::FieldVector fieldVector = stk::io::get_transient_fields(mBrokerB.meta_data());

  aTransientFields.assign(fieldVector.begin(), fieldVector.end());

  mBrokerA.get_global_variable_names(aGlobalVariableNames);

  for (const std::string& globalVariableName : aGlobalVariableNames) {
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
  for (TransientTransferByIdForRank *transfer : mTransfers) {
    transfer->do_transfer();
  }
}

size_t TransientFieldTransferById::setup_output_transient_fields(const std::string &parallelOutputMeshName)
{
  size_t outputFileIndex = mBrokerB.create_output_mesh(parallelOutputMeshName, stk::io::WRITE_RESULTS);
  mBrokerB.set_reference_input_region(outputFileIndex, mBrokerA);

  stk::mesh::FieldVector allTransientFieldsB = stk::io::get_transient_fields(mBrokerB.meta_data());

  for (stk::mesh::FieldBase* field : allTransientFieldsB) {
    mBrokerB.add_field(outputFileIndex, *field);
  }

  return outputFileIndex;
}

void TransientFieldTransferById::initialize(const std::vector<stk::mesh::EntityRank>& entityRanks)
{
  for (stk::mesh::EntityRank rank : entityRanks) {
    if (stk::io::get_transient_fields(mBrokerA.meta_data(), rank).size() > 0) {
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

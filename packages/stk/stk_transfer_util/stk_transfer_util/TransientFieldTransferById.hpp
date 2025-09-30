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

#include "stk_transfer/copy_by_id/TransferCopyById.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdStkMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/SearchByIdGeometric.hpp"

#include <stk_io/StkMeshIoBroker.hpp>

namespace Ioss { class Region; }

namespace stk {
namespace mesh { class BulkData; }
namespace transfer_utils {

class RepeatedTransferCopyByIdStkMeshAdapter : public stk::transfer::TransferCopyByIdStkMeshAdapter
{
public:
  RepeatedTransferCopyByIdStkMeshAdapter(stk::mesh::BulkData & mesh,
                                         const EntityVector & entities,
                                         const FieldVector & fields)
    : RepeatedTransferCopyByIdStkMeshAdapter(mesh, entities, fields, mesh.parallel())
  {}

  RepeatedTransferCopyByIdStkMeshAdapter(stk::mesh::BulkData & mesh,
                                         const EntityVector & entities,
                                         const FieldVector & fields,
                                         stk::ParallelMachine global_comm)
    : TransferCopyByIdStkMeshAdapter(mesh, entities, fields, global_comm),
      m_areValuesCached(false)
  {
    m_cachedFieldValues.resize(fields.size());
    m_cachedFieldValueIndex.resize(fields.size(), 0);
    m_cachedFieldValueSizes.resize(fields.size());
    m_cachedFieldValueSizeIndex.resize(fields.size(), 0);
  }

  virtual ~RepeatedTransferCopyByIdStkMeshAdapter() override = default;

  virtual void begin_transfer() const override
  {
    if (!m_areValuesCached) {
      for (auto & fieldValues : m_cachedFieldValues) {
        fieldValues.clear();
      }
      for (auto & fieldValueSizes : m_cachedFieldValueSizes) {
        fieldValueSizes.clear();
      }
    }
    std::fill(m_cachedFieldValueIndex.begin(), m_cachedFieldValueIndex.end(), 0);
    std::fill(m_cachedFieldValueSizeIndex.begin(), m_cachedFieldValueSizeIndex.end(), 0);
  }

  virtual void end_transfer() const override
  {
    m_areValuesCached = true;
  }

  void* get_cached_field_data_value(const stk::mesh::EntityKey & key, unsigned field_index) const
  {
    if (m_areValuesCached) {
      return m_cachedFieldValues[field_index][m_cachedFieldValueIndex[field_index]++];
    }
    else {
      const mesh::Entity entity = m_mesh.get_entity(key);
      stk::mesh::FieldBase* field = m_transfer_fields[field_index];
      void* fieldData = stk::mesh::field_data(*field, entity);
      m_cachedFieldValues[field_index].push_back(fieldData);
      return fieldData;
    }
  }

  virtual const void* field_data(const Mesh_ID & id, const unsigned field_index) const override
  {
    return get_cached_field_data_value(static_cast<EntityKey::entity_key_t>(id), field_index);
  }

  virtual void* field_data(const Mesh_ID & id, const unsigned field_index) override
  {
    return get_cached_field_data_value(static_cast<EntityKey::entity_key_t>(id), field_index);
  }

  unsigned get_cached_field_data_size(const stk::mesh::EntityKey & key, unsigned field_index) const
  {
    if (m_areValuesCached) {
      return m_cachedFieldValueSizes[field_index][m_cachedFieldValueSizeIndex[field_index]++];
    }
    else {
      const mesh::Entity entity = m_mesh.get_entity(key);
      stk::mesh::FieldBase* field = m_transfer_fields[field_index];
      unsigned fieldValueSize = stk::mesh::field_bytes_per_entity(*field, entity);
      m_cachedFieldValueSizes[field_index].push_back(fieldValueSize);
      return fieldValueSize;
    }
  }

  virtual unsigned field_data_size(const Mesh_ID & id, const unsigned field_index) const override
  {
    return get_cached_field_data_size(static_cast<stk::mesh::EntityKey::entity_key_t>(id), field_index);
  }

private:
  mutable bool m_areValuesCached;
  mutable std::vector<std::vector<void*>> m_cachedFieldValues;
  mutable std::vector<unsigned> m_cachedFieldValueIndex;
  mutable std::vector<std::vector<unsigned short>> m_cachedFieldValueSizes;
  mutable std::vector<unsigned> m_cachedFieldValueSizeIndex;
};

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

  RepeatedTransferCopyByIdStkMeshAdapter *mTransferMeshA = nullptr;
  RepeatedTransferCopyByIdStkMeshAdapter *mTransferMeshB = nullptr;

  stk::transfer::SearchByIdGeometric  mSearch;
  stk::transfer::TransferCopyById    *mTransfer = nullptr;

private:
  TransientTransferByIdForRank();

private:
  RepeatedTransferCopyByIdStkMeshAdapter *create_transfer_mesh(stk::mesh::MetaData &meta);
};

class TransientFieldTransferById
{
public:
  TransientFieldTransferById(stk::io::StkMeshIoBroker &brokerA, stk::io::StkMeshIoBroker &brokerB);

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

}
}

#endif // TRANSIENT_FIELD_TRANSFER_BY_ID_

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

#ifndef NGPFIELDPARALLEL_HPP
#define NGPFIELDPARALLEL_HPP

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_ngp/Ngp.hpp"
#include "stk_ngp/NgpField.hpp"
#include "stk_ngp/NgpParallelComm.hpp"
#include <vector>

namespace ngp {

template <typename T>
void parallel_sum(const stk::mesh::BulkData & bulk, const std::vector<ngp::Field<T> *> & ngpFields,
                  bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (ngp::Field<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::parallel_sum(bulk, stkFields);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void copy_owned_to_shared(const stk::mesh::BulkData & bulk, const std::vector<ngp::Field<T> *> & ngpFields,
                          bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (ngp::Field<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::copy_owned_to_shared(bulk, stkFields);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::Ghosting & ghosting,
                            const std::vector<ngp::Field<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = ghosting.mesh().mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (ngp::Field<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::communicate_field_data(ghosting, stkFields);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::BulkData & bulk,
                            const std::vector<ngp::Field<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (ngp::Field<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::communicate_field_data(bulk, stkFields);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
struct NgpFieldInfo
{
  NgpFieldInfo(ngp::Field<T>& fld)
    : m_field(fld) {}

  STK_FUNCTION
  NgpFieldInfo() = default;

  STK_FUNCTION
  NgpFieldInfo(const NgpFieldInfo&) = default;

  STK_FUNCTION
  ~NgpFieldInfo() = default;

  STK_FUNCTION
  operator const ngp::Field<T>&() const { return m_field; }

  ngp::Field<T> m_field;
};

template <typename T>
using FieldDataViewType = Kokkos::View<T*, MemSpace>;

template <typename T>
using FieldView = Kokkos::View<NgpFieldInfo<T>*, MemSpace>;

template <typename T>
class ParallelSumDataExchangeSymPackUnpackHandler
{
public:
  ParallelSumDataExchangeSymPackUnpackHandler(const ngp::Mesh& mesh, const std::vector<ngp::Field<T> *> & ngpFields)
    : m_ngpMesh(mesh),
      m_ngpFields(ngpFields),
      m_ngpFieldsOnDevice(FieldView<T>("ngpFieldsOnDevice", ngpFields.size()))
  {
    typename FieldView<T>::HostMirror ngpFieldsHostMirror = Kokkos::create_mirror_view(m_ngpFieldsOnDevice);
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFields.size(); fieldIdx++)
    {
      ngpFieldsHostMirror(fieldIdx) = NgpFieldInfo<T>(*m_ngpFields[fieldIdx]);
    }
    Kokkos::deep_copy(m_ngpFieldsOnDevice, ngpFieldsHostMirror);
  }

  ParallelSumDataExchangeSymPackUnpackHandler(const ParallelSumDataExchangeSymPackUnpackHandler & rhs) = default;

  void hostSizeMessages(int proc, size_t & numValues) const
  {
    numValues = 0;
    for (ngp::Field<T>* field : m_ngpFields)
    {
      stk::mesh::FieldBase* stkField = m_ngpMesh.get_bulk_on_host().mesh_meta_data().get_fields()[field->get_ordinal()];
      const stk::mesh::BucketIndices & stkBktIndices = m_ngpMesh.get_bulk_on_host().volatile_fast_shared_comm_map(field->get_rank())[proc];
      for (size_t i = 0; i < stkBktIndices.bucket_info.size(); ++i) {
        const unsigned bucketId = stkBktIndices.bucket_info[i].bucket_id;
        const unsigned numEntitiesThisBucket = stkBktIndices.bucket_info[i].num_entities_this_bucket;
        const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(*stkField, bucketId);
        numValues += numScalarsPerEntity * numEntitiesThisBucket;
      }
    }
  }

  STK_FUNCTION
  void devicePackMessage(int /*myProc*/, int proc, FieldDataViewType<T> & sendData) const
  {
    size_t fieldCount = 0;
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFieldsOnDevice.size(); ++fieldIdx)
    {
      const ngp::Field<T> & field = m_ngpFieldsOnDevice(fieldIdx);
      ngp::DeviceCommMapIndices deviceCommMapIndices = m_ngpMesh.volatile_fast_shared_comm_map(field.get_rank(), proc);
      for (size_t entry = 0; entry < deviceCommMapIndices.size(); ++entry) {
        size_t numComponents = field.get_num_components_per_entity(deviceCommMapIndices[entry]);
        for (size_t comp = 0; comp < numComponents; ++comp) {
          sendData(fieldCount++) = field.get(deviceCommMapIndices[entry], comp);
        }
      }
    }
  }

  STK_FUNCTION
  void deviceUnpackMessage(int /*myProc*/, int proc, FieldDataViewType<T> & recvData) const
  {
    size_t fieldCount = 0;
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFieldsOnDevice.size(); ++fieldIdx)
    {
      const ngp::Field<T> & field = m_ngpFieldsOnDevice(fieldIdx);
      ngp::DeviceCommMapIndices deviceCommMapIndices = m_ngpMesh.volatile_fast_shared_comm_map(field.get_rank(), proc);
      for (size_t entry = 0; entry < deviceCommMapIndices.size(); ++entry) {
        size_t numComponents = field.get_num_components_per_entity(deviceCommMapIndices[entry]);
        for (size_t comp = 0; comp < numComponents; ++comp) {
          field.get(deviceCommMapIndices[entry], comp) += recvData(fieldCount++);
        }
      }
    }
  }

private:
  const ngp::Mesh m_ngpMesh;
  const std::vector<ngp::Field<T> *>& m_ngpFields;
  FieldView<T> m_ngpFieldsOnDevice;
};

template <typename T>
void parallel_sum_device_mpi(const ngp::Mesh& ngpMesh, const std::vector<ngp::Field<T> *> & ngpFields)
{
  const stk::mesh::BulkData& bulk = ngpMesh.get_bulk_on_host();

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->sync_to_device();
  }

  ParallelSumDataExchangeSymPackUnpackHandler<T> exchangeHandler(ngpMesh, ngpFields);

  stk::mesh::EntityRank first_field_rank = ngpFields[0]->get_rank();
  std::vector<int> comm_procs = bulk.all_sharing_procs(first_field_rank);
  for (size_t i = 1; i < ngpFields.size(); ++i) {
    const ngp::Field<T> * f = ngpFields[i];
    if (f->get_rank() != first_field_rank) {
      const std::vector<int>& sharing_procs = bulk.all_sharing_procs(f->get_rank());
      for(int p : sharing_procs) {
        stk::util::insert_keep_sorted_and_unique(p, comm_procs);
      }
    }
  }

  const bool deterministic = false;
  ngp::parallel_data_exchange_sym_pack_unpack<double>(MPI_COMM_WORLD,
                                                      comm_procs,
                                                      exchangeHandler,
                                                      deterministic);

  for (ngp::Field<T> * ngpField : ngpFields) {
    ngpField->modify_on_device();
  }
}

}

#endif


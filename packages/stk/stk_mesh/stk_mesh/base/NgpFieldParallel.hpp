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

#ifndef STK_MESH_NGPFIELDPARALLEL_HPP
#define STK_MESH_NGPFIELDPARALLEL_HPP

#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Ngp.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/NgpParallelComm.hpp"
#include <vector>

namespace stk {
namespace mesh {

template <typename T>
void parallel_sum(const stk::mesh::BulkData & bulk, const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                  bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
    ngpField->sync_to_host();
  }

  stk::mesh::parallel_sum(bulk, stkFields);

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->modify_on_host();

    if (doFinalSyncBackToDevice) {
      ngpField->sync_to_device();
    }
  }
}

template <typename T>
void do_final_sync_to_device(const std::vector<NgpField<T>*>& ngpFields)
{
  for(auto ngpField : ngpFields) {
    ngpField->sync_to_device();
  }
}

template <typename T>
void copy_owned_to_shared(const stk::mesh::BulkData & bulk,
                          const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                          bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::copy_owned_to_shared(bulk, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::Ghosting & ghosting,
                            const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = ghosting.mesh().mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::communicate_field_data(ghosting, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
void communicate_field_data(const stk::mesh::BulkData & bulk,
                            const std::vector<stk::mesh::NgpField<T> *> & ngpFields,
                            bool doFinalSyncBackToDevice = true)
{
  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  const std::vector<stk::mesh::FieldBase *> & allStkFields = meta.get_fields();

  std::vector<const stk::mesh::FieldBase *> stkFields;
  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    stkFields.push_back(allStkFields[ngpField->get_ordinal()]);
  }

  stk::mesh::communicate_field_data(bulk, stkFields);

  if(doFinalSyncBackToDevice) {
    do_final_sync_to_device(ngpFields);
  }
}

template <typename T>
struct NgpFieldInfo
{
  NgpFieldInfo(stk::mesh::NgpField<T>& fld)
    : m_field(fld) {}

  KOKKOS_DEFAULTED_FUNCTION
  NgpFieldInfo() = default;

  KOKKOS_DEFAULTED_FUNCTION
  NgpFieldInfo(const NgpFieldInfo&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~NgpFieldInfo() = default;

  KOKKOS_FUNCTION
  operator const stk::mesh::NgpField<T>&() const { return m_field; }

  stk::mesh::NgpField<T> m_field;
};

template <typename T>
using FieldDataViewType = Kokkos::View<T*, stk::ngp::MemSpace>;

template <typename T>
using FieldView = Kokkos::View<NgpFieldInfo<T>*, stk::ngp::MemSpace>;

template <typename T>
class ParallelSumDataExchangeSymPackUnpackHandler
{
public:
  ParallelSumDataExchangeSymPackUnpackHandler(const stk::mesh::NgpMesh& mesh, const std::vector<stk::mesh::NgpField<T> *> & ngpFields)
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
    for (stk::mesh::NgpField<T>* field : m_ngpFields)
    {
      stk::mesh::FieldBase* stkField = m_ngpMesh.get_bulk_on_host().mesh_meta_data().get_fields()[field->get_ordinal()];
      stk::mesh::HostCommMapIndices  commMapIndices = m_ngpMesh.get_bulk_on_host().volatile_fast_shared_comm_map(field->get_rank(), proc);
      for (size_t i = 0; i < commMapIndices.extent(0); ++i) {
        const unsigned bucketId = commMapIndices(i).bucket_id;
        const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(*stkField, bucketId);
        numValues += numScalarsPerEntity;
      }
    }
  }

  KOKKOS_FUNCTION
  void devicePackMessage(int /*myProc*/, int proc, FieldDataViewType<T> & sendData) const
  {
    size_t fieldCount = 0;
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFieldsOnDevice.size(); ++fieldIdx)
    {
      const stk::mesh::NgpField<T> & field = m_ngpFieldsOnDevice(fieldIdx);
      stk::mesh::DeviceCommMapIndices deviceCommMapIndices = m_ngpMesh.volatile_fast_shared_comm_map(field.get_rank(), proc);
      for (size_t entry = 0; entry < deviceCommMapIndices.size(); ++entry) {
        size_t numComponents = field.get_num_components_per_entity(deviceCommMapIndices[entry]);
        for (size_t comp = 0; comp < numComponents; ++comp) {
          sendData(fieldCount++) = field.get(deviceCommMapIndices[entry], comp);
        }
      }
    }
  }

  KOKKOS_FUNCTION
  void deviceUnpackMessage(int /*myProc*/, int proc, FieldDataViewType<T> & recvData) const
  {
    size_t fieldCount = 0;
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFieldsOnDevice.size(); ++fieldIdx)
    {
      const stk::mesh::NgpField<T> & field = m_ngpFieldsOnDevice(fieldIdx);
      stk::mesh::DeviceCommMapIndices deviceCommMapIndices = m_ngpMesh.volatile_fast_shared_comm_map(field.get_rank(), proc);
      for (size_t entry = 0; entry < deviceCommMapIndices.size(); ++entry) {
        size_t numComponents = field.get_num_components_per_entity(deviceCommMapIndices[entry]);
        for (size_t comp = 0; comp < numComponents; ++comp) {
          field.get(deviceCommMapIndices[entry], comp) += recvData(fieldCount++);
        }
      }
    }
  }

private:
  const stk::mesh::NgpMesh m_ngpMesh;
  const std::vector<stk::mesh::NgpField<T> *>& m_ngpFields;
  FieldView<T> m_ngpFieldsOnDevice;
};

template <typename T>
void parallel_sum_device_mpi(const stk::mesh::NgpMesh& ngpMesh, const std::vector<stk::mesh::NgpField<T> *> & ngpFields)
{
  const stk::mesh::BulkData& bulk = ngpMesh.get_bulk_on_host();

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->sync_to_device();
  }

  ParallelSumDataExchangeSymPackUnpackHandler<T> exchangeHandler(ngpMesh, ngpFields);

  stk::mesh::EntityRank first_field_rank = ngpFields[0]->get_rank();
  std::vector<int> comm_procs = bulk.all_sharing_procs(first_field_rank);
  for (size_t i = 1; i < ngpFields.size(); ++i) {
    const stk::mesh::NgpField<T> * f = ngpFields[i];
    if (f->get_rank() != first_field_rank) {
      const std::vector<int>& sharing_procs = bulk.all_sharing_procs(f->get_rank());
      for(int p : sharing_procs) {
        stk::util::insert_keep_sorted_and_unique(p, comm_procs);
      }
    }
  }

  const bool deterministic = false;
  stk::mesh::ngp_parallel_data_exchange_sym_pack_unpack<double>(bulk.parallel(),
                                                                comm_procs,
                                                                exchangeHandler,
                                                                deterministic);

  for (stk::mesh::NgpField<T> * ngpField : ngpFields) {
    ngpField->modify_on_device();
  }
}

}
}

#endif


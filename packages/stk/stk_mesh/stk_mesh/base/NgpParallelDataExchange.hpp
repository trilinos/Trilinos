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

#ifndef STK_MESH_NGP_PARALLEL_DATA_EXCHANGE_HPP
#define STK_MESH_NGP_PARALLEL_DATA_EXCHANGE_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/MPI.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {

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
using BufferViewType = Kokkos::View<T*, stk::ngp::MemSpace>;

template <typename T>
class ParallelSumDataExchangeSymPackUnpackHandler
{
public:
  ParallelSumDataExchangeSymPackUnpackHandler(const stk::mesh::NgpMesh& mesh, const std::vector<stk::mesh::NgpField<T> *> & ngpFields)
    : m_ngpMesh(const_cast<stk::mesh::NgpMesh&>(mesh)),
      m_ngpFields(ngpFields),
      m_ngpFieldsOnDevice(FieldView<T>("ngpFieldsOnDevice", ngpFields.size())),
      m_deviceSendData(BufferViewType<T>("deviceSendData", 1)),
      m_deviceRecvData(BufferViewType<T>("deviceRecvData", 1))
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
      stk::mesh::HostCommMapIndices  commMapIndices = m_ngpMesh.get_bulk_on_host().template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(field->get_rank(), proc);
      for (size_t i = 0; i < commMapIndices.extent(0); ++i) {
        const unsigned bucketId = commMapIndices(i).bucket_id;
        const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(*stkField, bucketId);
        numValues += numScalarsPerEntity;
      }
    }
  }

  stk::mesh::NgpMesh& get_ngp_mesh() const {
    return m_ngpMesh;
  }

  std::vector<stk::mesh::NgpField<T> *> const& get_ngp_fields() const {
    return m_ngpFields;
  }

  KOKKOS_FUNCTION
  FieldView<T> get_ngp_fields_on_device() const {
    return m_ngpFieldsOnDevice;
  }

  BufferViewType<T> get_device_send_data() const {
    return m_deviceSendData;
  }

  BufferViewType<T> get_device_recv_data() const {
    return m_deviceRecvData;
  }

  UnsignedViewType::HostMirror get_host_buffer_offsets() const {
    return m_ngpMesh.get_ngp_parallel_sum_host_buffer_offsets();
  }

  Unsigned2dViewType::HostMirror get_host_mesh_indices_offsets() const {
    return m_ngpMesh.get_ngp_parallel_sum_host_mesh_indices_offsets();
  }

  Unsigned2dViewType get_device_mesh_indices_offsets() const {
    return m_ngpMesh.get_ngp_parallel_sum_device_mesh_indices_offsets();
  }

  void resize_device_mpi_buffers(size_t size) {
    Kokkos::resize(Kokkos::WithoutInitializing, m_deviceSendData, size);
    Kokkos::resize(Kokkos::WithoutInitializing, m_deviceRecvData, size);
  }

private:
  stk::mesh::NgpMesh& m_ngpMesh;
  const std::vector<stk::mesh::NgpField<T> *>& m_ngpFields;
  FieldView<T> m_ngpFieldsOnDevice;

  BufferViewType<T> m_deviceSendData;
  BufferViewType<T> m_deviceRecvData;
};

}
}

#endif

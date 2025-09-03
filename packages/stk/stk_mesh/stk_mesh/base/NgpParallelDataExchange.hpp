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
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {

namespace impl {

template<typename NgpFieldType>
void fill_field_ranks(const std::vector<NgpFieldType*>& fields,
                      std::vector<stk::mesh::EntityRank>& fieldRanks)
{
  fieldRanks.clear();
  fieldRanks.reserve(fields.size());
  for(const NgpFieldType* f : fields) {
    fieldRanks.push_back(f->get_rank());
  }
  stk::util::sort_and_unique(fieldRanks);
}

} // namespace impl

template <typename NgpFieldType>
struct NgpFieldInfo
{
  NgpFieldInfo(NgpFieldType& fld)
    : m_field(fld) {}

  KOKKOS_DEFAULTED_FUNCTION
  NgpFieldInfo() = default;

  KOKKOS_DEFAULTED_FUNCTION
  NgpFieldInfo(const NgpFieldInfo&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~NgpFieldInfo() = default;

  KOKKOS_FUNCTION
  operator const NgpFieldType&() const { return m_field; }

  NgpFieldType m_field;
};

template <typename T, typename MEMSPACE>
using FieldDataViewType = Kokkos::View<T*, MEMSPACE>;

template <typename NgpFieldType, typename MEMSPACE>
using FieldView = Kokkos::View<NgpFieldInfo<NgpFieldType>*, MEMSPACE>;

template <typename T, typename MEMSPACE>
using BufferViewType = Kokkos::View<T*, MEMSPACE>;

template <typename NgpMeshType, typename NgpFieldType>
class ParallelSumDataExchangeSymPackUnpackHandler
{
public:
  using T = typename NgpFieldType::value_type;
  using mem_space = typename NgpMeshType::ngp_mem_space;
  using mesh_type = NgpMeshType;
  using field_type = NgpFieldType;

  ParallelSumDataExchangeSymPackUnpackHandler(const NgpMeshType& mesh, const std::vector<NgpFieldType *> & ngpFields)
    : m_ngpMesh(const_cast<NgpMeshType&>(mesh)),
      m_fieldRanks(),
      m_ngpFields(ngpFields),
      m_ngpFieldsOnDevice(FieldView<NgpFieldType,mem_space>("ngpFieldsOnDevice", ngpFields.size())),
      m_deviceSendData(BufferViewType<T,mem_space>("deviceSendData", 1)),
      m_deviceRecvData(BufferViewType<T,mem_space>("deviceRecvData", 1))
  {
    impl::fill_field_ranks(ngpFields, m_fieldRanks);
    typename FieldView<NgpFieldType,mem_space>::host_mirror_type ngpFieldsHostMirror = Kokkos::create_mirror_view(m_ngpFieldsOnDevice);
    for (size_t fieldIdx = 0; fieldIdx < m_ngpFields.size(); fieldIdx++)
    {
      ngpFieldsHostMirror(fieldIdx) = NgpFieldInfo<NgpFieldType>(*m_ngpFields[fieldIdx]);
    }
    Kokkos::deep_copy(m_ngpFieldsOnDevice, ngpFieldsHostMirror);
  }

  ParallelSumDataExchangeSymPackUnpackHandler(const ParallelSumDataExchangeSymPackUnpackHandler & rhs) = default;

  void hostSizeMessages(int proc, size_t & numValues, bool includeGhosts=false) const
  {
    numValues = 0;
    for(stk::mesh::EntityRank fieldRank : m_fieldRanks) {
      stk::mesh::HostCommMapIndices<stk::ngp::MemSpace> commMapIndices =
        m_ngpMesh.get_bulk_on_host().template volatile_fast_shared_comm_map<stk::ngp::MemSpace>(fieldRank, proc, includeGhosts);

      for (size_t i = 0; i < commMapIndices.extent(0); ++i) {
        const unsigned bucketId = commMapIndices(i).bucket_id;

        for (NgpFieldType* field : m_ngpFields) {
          if (field->get_rank() == fieldRank) {
            stk::mesh::FieldBase* stkField = m_ngpMesh.get_bulk_on_host().mesh_meta_data().get_fields()[field->get_ordinal()];
            const unsigned numScalarsPerEntity = stk::mesh::field_scalars_per_entity(*stkField, bucketId);
            numValues += numScalarsPerEntity;
          }
        }
      }
    }
  }

  NgpMeshType& get_ngp_mesh() const {
    return m_ngpMesh;
  }

  const std::vector<stk::mesh::EntityRank>& get_field_ranks() const {
    return m_fieldRanks;
  }

  std::vector<NgpFieldType*> const& get_ngp_fields() const {
    return m_ngpFields;
  }

  KOKKOS_FUNCTION
  const FieldView<NgpFieldType,mem_space>& get_ngp_fields_on_device() const
  {
    return m_ngpFieldsOnDevice;
  }

  BufferViewType<T,mem_space>& get_device_send_data() {
    return m_deviceSendData;
  }

  BufferViewType<T,mem_space>& get_device_recv_data() {
    return m_deviceRecvData;
  }

  typename UnsignedViewType<stk::ngp::MemSpace>::host_mirror_type& get_host_buffer_offsets()
  {
    return m_ngpMesh.get_ngp_parallel_sum_host_buffer_offsets();
  }

  typename UnsignedViewType<stk::ngp::MemSpace>::host_mirror_type& get_host_mesh_indices_offsets()
  {
    return m_ngpMesh.get_ngp_parallel_sum_host_mesh_indices_offsets();
  }

  auto& get_device_mesh_indices_offsets()
  {
    return m_ngpMesh.get_ngp_parallel_sum_device_mesh_indices_offsets();
  }

  void resize_device_mpi_buffers(size_t size) {
    if (m_deviceSendData.size() < size) {
      Kokkos::resize(Kokkos::WithoutInitializing, m_deviceSendData, size);
    }
    if (m_deviceRecvData.size() < size) {
      Kokkos::resize(Kokkos::WithoutInitializing, m_deviceRecvData, size);
    }
  }

private:
  NgpMeshType& m_ngpMesh;
  std::vector<stk::mesh::EntityRank> m_fieldRanks;
  const std::vector<NgpFieldType *>& m_ngpFields;
  FieldView<NgpFieldType,mem_space> m_ngpFieldsOnDevice;

  BufferViewType<T,mem_space> m_deviceSendData;
  BufferViewType<T,mem_space> m_deviceRecvData;
};

}
}

#endif

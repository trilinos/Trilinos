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

#ifndef STK_MESH_NGPUTILS_HPP
#define STK_MESH_NGPUTILS_HPP

#include <stk_util/stk_config.h>
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Types.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_util/util/StkNgpVector.hpp>
#include <numeric>

namespace stk {
namespace mesh {

inline void ngp_field_fence(MetaData& meta)
{
  auto fields = meta.get_fields();

  for(auto field : fields) {
    if (field->has_device_data()) {
      field->fence();
    }
  }
}

inline void require_ngp_mesh_rank_limit(const stk::mesh::MetaData& meta)
{
  const size_t maxNumRanks = stk::topology::NUM_RANKS;
  const size_t numRanks = meta.entity_rank_count();
  STK_ThrowRequireMsg(numRanks <= maxNumRanks,
                "stk::mesh::NgpMesh: too many entity ranks ("<<numRanks
                <<"). Required to be less-or-equal stk::topology::NUM_RANKS");
}

inline stk::NgpVector<unsigned> get_bucket_ids(const stk::mesh::BulkData &bulk,
                                               stk::mesh::EntityRank rank,
                                               const stk::mesh::Selector &selector)
{
  Kokkos::Profiling::pushRegion("get_bucket_ids");
  const stk::mesh::BucketVector &buckets = bulk.get_buckets(rank, selector);
  stk::NgpVector<unsigned> bucketIds(buckets.size());
  for(size_t i=0; i<buckets.size(); i++) {
    bucketIds[i] = buckets[i]->bucket_id();
  }
  bucketIds.copy_host_to_device();
  Kokkos::Profiling::popRegion();
  return bucketIds;
}

template <typename NgpSpace, typename FieldDataBytesViewType>
void assemble_field_data_bytes_on_device(const std::vector<FieldBase*>& fields, FieldDataBytesViewType& fieldDataBytesView)
{
  Kokkos::resize(Kokkos::WithoutInitializing, fieldDataBytesView, fields.size());
  auto fieldDataBytesViewOnHost = Kokkos::create_mirror(fieldDataBytesView);
  for (size_t fieldIdx = 0; fieldIdx < fields.size(); ++fieldIdx) {
    const auto& field = *fields[fieldIdx];
    fieldDataBytesViewOnHost(fieldIdx) = field.data_bytes<std::byte, NgpSpace>();
  }
  Kokkos::deep_copy(fieldDataBytesView, fieldDataBytesViewOnHost);
}

struct BytePtr
{
  std::byte* ptr;
};

template<typename InitValsPtrViewType>
inline
void assemble_field_init_vals_on_device(const std::vector<FieldBase*>& fields, InitValsPtrViewType& initValsPtrsView)
{
  Kokkos::resize(initValsPtrsView, fields.size());
  for (unsigned fieldIdx = 0; fieldIdx < fields.size(); ++fieldIdx) {
    const auto& initVals = fields[fieldIdx]->get_initial_value_bytes();
    initValsPtrsView(fieldIdx) = BytePtr{initVals.extent(0) > 0 ? initVals.data() : nullptr};
  }
}

template<stk::mesh::Layout LayoutValue,
         typename SrcFieldBytesType, typename DestFieldBytesType>
KOKKOS_INLINE_FUNCTION
void copy_bytes_kernel(const SrcFieldBytesType& srcFieldBytes,
                             DestFieldBytesType& destFieldBytes,
                       const FastMeshIndex& srcIndex,
                       const FastMeshIndex& destIndex,
                       const std::byte* initVals)
{
  auto srcBytes = srcFieldBytes.template entity_bytes<LayoutValue>(srcIndex);
  auto destBytes = destFieldBytes.template entity_bytes<LayoutValue>(destIndex);
  if (srcBytes.num_bytes() > 0) {
    for(ByteIdx idx : destBytes.bytes()) {
      destBytes(idx) = srcBytes(idx);
    }
  }
  else {
    if (initVals != nullptr) {
      for(ByteIdx idx : destBytes.bytes()) {
        destBytes(idx) = initVals[idx()];
      }
    }
    else {
      for(ByteIdx idx : destBytes.bytes()) {
        destBytes(idx) = static_cast<std::byte>(0);
      }
    }
  }
}

template<stk::mesh::Layout LayoutValue,
         typename DeviceFieldDataBytesViewType, typename InitValsViewType>
KOKKOS_INLINE_FUNCTION
void copy_entity_bytes_kernel(DeviceFieldDataBytesViewType& deviceFieldDataBytes,
                              const FastMeshIndex& srcIndex,
                              const FastMeshIndex& destIndex,
                              const InitValsViewType& initValsView)
{
  for(unsigned i=0; i<deviceFieldDataBytes.extent(0); ++i) {
    copy_bytes_kernel<Layout::Left>(deviceFieldDataBytes(i), deviceFieldDataBytes(i), srcIndex, destIndex, initValsView(i).ptr);
  }
}

template<typename Space>
void copy_entity_field_bytes(const std::vector<FieldBase*>& fields,
                        const FastMeshIndex& srcIndex, const FastMeshIndex& destIndex)
{
  STK_ThrowErrorMsg("un-specialized copy_entity_field_bytes template, call with either stk::ngp::HostSpace or stk::ngp::DeviceSpace");
}

template<>
inline
void copy_entity_field_bytes<stk::ngp::HostSpace>(const std::vector<FieldBase*>& fields,
                        const FastMeshIndex& srcIndex, const FastMeshIndex& destIndex)
{
  for(FieldBase* field : fields) {
    auto& srcFieldBytes = field->data_bytes<const std::byte,stk::ngp::HostSpace>();
    auto& destFieldBytes = field->data_bytes<std::byte,stk::ngp::HostSpace>();
    const auto& initVals = field->get_initial_value_bytes();
    const std::byte* initValsPtr = initVals.extent(0) > 0 ? initVals.data() : nullptr;
    if (field->host_data_layout() == Layout::Left) {
      copy_bytes_kernel<Layout::Left>(srcFieldBytes, destFieldBytes, srcIndex, destIndex, initValsPtr);
    }
    else if (field->host_data_layout() == Layout::Right) {
      copy_bytes_kernel<Layout::Right>(srcFieldBytes, destFieldBytes, srcIndex, destIndex, initValsPtr);
    }
    else {
      STK_ThrowErrorMsg("Unsupported host Field data layout: " << field->host_data_layout());
    }
  }
}

#ifdef STK_USE_DEVICE_MESH

template<>
inline
void copy_entity_field_bytes<stk::ngp::DeviceSpace>(const std::vector<FieldBase*>& fields,
                        const FastMeshIndex& srcIndex, const FastMeshIndex& destIndex)
{
  auto notParallel = Kokkos::RangePolicy<stk::ngp::DeviceSpace::exec_space>(0, 1);

  for(FieldBase* field : fields) {
    auto& srcFieldBytes = field->data_bytes<const std::byte,stk::ngp::DeviceSpace>();
    auto& destFieldBytes = field->data_bytes<std::byte,stk::ngp::DeviceSpace>();
    const auto& initVals = field->get_initial_value_bytes();
    const std::byte* initValsPtr = initVals.extent(0) > 0 ? initVals.data() : nullptr;
    STK_ThrowAssertMsg(field->device_data_layout() == Layout::Left,"Device field layout must always be Left");
    Kokkos::parallel_for("copy_entity_fields Left", notParallel,
      KOKKOS_LAMBDA(const int) {
        copy_bytes_kernel<Layout::Left>(srcFieldBytes, destFieldBytes, srcIndex, destIndex, initValsPtr);
      }
    );
  }
}

#endif // STK_USE_DEVICE_MESH

}
}

#endif


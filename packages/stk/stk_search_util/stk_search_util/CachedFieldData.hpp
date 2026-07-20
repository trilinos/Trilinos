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
//

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDFIELDDATA_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDFIELDDATA_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/FieldData.hpp"
#include "stk_mesh/base/FieldDataBytes.hpp"
#include "stk_mesh/base/FieldState.hpp"  // for FieldState, StateNP1, StateNM1
#include "stk_mesh/base/Types.hpp"       // for EntityRank
#include "stk_util/util/string_case_compare.hpp"

#include <functional>                    // for function
#include <limits>                        // for numeric_limits
#include <string>                        // for string, operator==, basic_st...
#include <utility>                       // for pair
#include <vector>                        // for vector
#include <type_traits>

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

class CachedEntityFieldDataBase
{
public:
  virtual ~CachedEntityFieldDataBase() = default;

  int componentStride{0};
  int numComponents{0};

  int copyStride{0};
  int numCopies{0};
};


template< typename DataType = double>
class CachedEntityFieldData : public CachedEntityFieldDataBase
{
public:
  DataType * pointer{nullptr};
  const DataType * constPointer{nullptr};
};

struct CachedFieldDataBase {
  CachedFieldDataBase() = default;
  virtual ~CachedFieldDataBase() = default;

  CachedFieldDataBase(const stk::mesh::FieldBase* field, unsigned fieldIndex)
  : m_field(field), m_fieldIndex(fieldIndex) {}

  virtual void populate_entity_data(stk::mesh::Entity /*entity*/, CachedEntityFieldDataBase& /*data*/) const = 0;

  virtual int component_stride(stk::mesh::Entity /*entity*/) const = 0;
  virtual int num_components(stk::mesh::Entity /*entity*/) const = 0;
  virtual int copy_stride(stk::mesh::Entity /*entity*/) const = 0;
  virtual int num_copies(stk::mesh::Entity /*entity*/) const = 0;

  virtual void* opaque_pointer(stk::mesh::Entity /*entity*/) = 0;
  virtual const void* opaque_pointer(stk::mesh::Entity /*entity*/) const = 0;

  virtual int bytes_per_scalar(stk::mesh::Entity /*entity*/) const = 0;
  virtual int scalar_byte_stride(stk::mesh::Entity /*entity*/) const = 0;
  virtual int num_bytes(stk::mesh::Entity /*entity*/) const = 0;

  const std::string& name() const {return m_field->name();}

  const stk::mesh::FieldBase* m_field{nullptr};
  unsigned m_fieldIndex{0};
};

template< typename DataType = double, stk::mesh::Layout DataLayout = stk::mesh::Layout::Right>
struct CachedFieldData : public CachedFieldDataBase {
  using BaseClass = CachedFieldDataBase;

  CachedFieldData() = default;
  ~CachedFieldData() = default;

  CachedFieldData(const stk::mesh::FieldData<DataType, stk::ngp::HostSpace, DataLayout>& fieldData,
                  const stk::mesh::FieldBase* field, unsigned fieldIndex)
  : CachedFieldDataBase(field, fieldIndex), m_fieldData(fieldData), m_fieldBytes(field->data_bytes<std::byte>()) {}

  void populate_entity_data(stk::mesh::Entity entity, CachedEntityFieldDataBase& data_) const override
  {
    CachedEntityFieldData<DataType>* data = dynamic_cast<CachedEntityFieldData<DataType>*>(&data_);
    STK_ThrowAssert(nullptr != data);

    auto fieldValues = m_fieldData.entity_values(entity);
    data->pointer = fieldValues.pointer();
    data->constPointer = data->pointer;
    data->componentStride = fieldValues.component_stride();
    data->numComponents = fieldValues.num_components();
    data->copyStride = fieldValues.copy_stride();
    data->numCopies = fieldValues.num_copies();
  }

  int component_stride(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.component_stride();
  }

  int num_components(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.num_components();
  }

  int copy_stride(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.copy_stride();
  }

  int num_copies(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.num_copies();
  }

  void* opaque_pointer(stk::mesh::Entity entity) override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return reinterpret_cast<void*>(fieldValues.pointer());
  }

  const void* opaque_pointer(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return reinterpret_cast<const void*>(fieldValues.pointer());
  }

  int bytes_per_scalar(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int bytesPerScalar = entityBytes.bytes_per_scalar();

    return bytesPerScalar;
  }

  int scalar_byte_stride(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int scalarByteStride = entityBytes.scalar_byte_stride();

    return scalarByteStride;
  }

  int num_bytes(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int numBytes = entityBytes.num_bytes();

    return numBytes;
  }

  stk::mesh::FieldData<DataType, stk::ngp::HostSpace, DataLayout> m_fieldData;
  stk::mesh::FieldDataBytes<stk::ngp::HostSpace> m_fieldBytes;
};

template< typename DataType = double,  stk::mesh::Layout DataLayout = stk::mesh::Layout::Right>
struct ConstCachedFieldData : public CachedFieldDataBase  {
  using BaseClass = CachedFieldDataBase;

  ConstCachedFieldData() = default;
  ~ConstCachedFieldData() = default;

  ConstCachedFieldData(const stk::mesh::ConstFieldData<DataType, stk::ngp::HostSpace, DataLayout>& fieldData,
                       const stk::mesh::FieldBase* field, unsigned fieldIndex)
  : CachedFieldDataBase(field, fieldIndex), m_fieldData(fieldData), m_fieldBytes(field->data_bytes<const std::byte>()) {}

  void populate_entity_data(stk::mesh::Entity entity, CachedEntityFieldDataBase& data_) const override
  {
    CachedEntityFieldData<DataType>* data = dynamic_cast<CachedEntityFieldData<DataType>*>(&data_);
    STK_ThrowAssert(nullptr != data);

    auto fieldValues = m_fieldData.entity_values(entity);
    data->pointer = nullptr;
    data->constPointer = fieldValues.pointer();
    data->componentStride = fieldValues.component_stride();
    data->numComponents = fieldValues.num_components();
    data->copyStride = fieldValues.copy_stride();
    data->numCopies = fieldValues.num_copies();
  }

  int component_stride(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.component_stride();
  }

  int num_components(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.num_components();
  }

  int copy_stride(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.copy_stride();
  }

  int num_copies(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return fieldValues.num_copies();
  }

  void* opaque_pointer(stk::mesh::Entity /*entity*/) override { return nullptr; }

  const void* opaque_pointer(stk::mesh::Entity entity) const override
  {
    auto fieldValues = m_fieldData.entity_values(entity);
    return reinterpret_cast<const void*>(fieldValues.pointer());
  }

  virtual int bytes_per_scalar(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int bytesPerScalar = entityBytes.bytes_per_scalar();

    return bytesPerScalar;
  }

  virtual int scalar_byte_stride(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int scalarByteStride = entityBytes.scalar_byte_stride();

    return scalarByteStride;
  }

  int num_bytes(stk::mesh::Entity entity) const override
  {
    auto entityBytes = m_fieldBytes.template entity_bytes<DataLayout>(entity);
    const int numBytes = entityBytes.num_bytes();

    return numBytes;
  }

  stk::mesh::ConstFieldData<DataType, stk::ngp::HostSpace, DataLayout> m_fieldData;
  stk::mesh::ConstFieldDataBytes<stk::ngp::HostSpace> m_fieldBytes;
};

std::shared_ptr<CachedFieldDataBase> get_cached_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex = 0);
void fill_cached_field_data(const stk::mesh::FieldBase* field, std::shared_ptr<CachedFieldDataBase>& cachedFieldData);
void fill_cached_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex, std::shared_ptr<CachedFieldDataBase>& cachedFieldData);
void fill_cached_field_data(const std::vector<stk::mesh::FieldBase*>& fieldVec, std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData);

std::shared_ptr<CachedFieldDataBase> get_cached_const_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex = 0);
void fill_cached_const_field_data(const stk::mesh::FieldBase* field, std::shared_ptr<CachedFieldDataBase>& cachedFieldData);
void fill_cached_const_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex, std::shared_ptr<CachedFieldDataBase>& cachedFieldData);
void fill_cached_const_field_data(const std::vector<stk::mesh::FieldBase*>& fieldVec, std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData);

void clear_cached_field_data(std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData);
void clear_cached_field_data(std::shared_ptr<CachedFieldDataBase>& cachedFieldData);

} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CACHEDFIELDDATA_HPP_ */

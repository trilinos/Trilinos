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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_search_util/CachedFieldData.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

namespace impl {
template<typename T>
void internal_fill_cached_field_data(const stk::mesh::FieldBase* field,
                                     const unsigned fieldIndex,
                                     std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  STK_ThrowAssert(nullptr != field);

  if(field->host_data_layout() == stk::mesh::Layout::Left) {
    auto fieldData = field->data<T, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>();
    cachedFieldData = std::make_shared< CachedFieldData<T, stk::mesh::Layout::Left> >(fieldData, field, fieldIndex);
  } else if(field->host_data_layout() == stk::mesh::Layout::Right) {
    auto fieldData = field->data<T, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>();
    cachedFieldData = std::make_shared< CachedFieldData<T, stk::mesh::Layout::Right> >(fieldData, field, fieldIndex);
  }
}

template<typename T>
void internal_fill_cached_const_field_data(const stk::mesh::FieldBase* field,
                                           const unsigned fieldIndex,
                                           std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  STK_ThrowAssert(nullptr != field);

  if(field->host_data_layout() == stk::mesh::Layout::Left) {
    auto fieldData = field->data<T, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>();
    cachedFieldData = std::make_shared< ConstCachedFieldData<T, stk::mesh::Layout::Left> >(fieldData, field, fieldIndex);
  } else if(field->host_data_layout() == stk::mesh::Layout::Right) {
    auto fieldData = field->data<T, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>();
    cachedFieldData = std::make_shared< ConstCachedFieldData<T, stk::mesh::Layout::Right> >(fieldData, field, fieldIndex);
  }
}
}

void fill_cached_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex, std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  if(nullptr == field) return;

  field_datatype_execute(*field,
    [&]<typename T>(const stk::mesh::FieldBase& /*fieldBase*/) {
      impl::internal_fill_cached_field_data<T>(field, fieldIndex, cachedFieldData);
    }
  );
}

void fill_cached_field_data(const stk::mesh::FieldBase* field, std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  fill_cached_field_data(field, 0u, cachedFieldData);
}

std::shared_ptr<CachedFieldDataBase> get_cached_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex)
{
  std::shared_ptr<CachedFieldDataBase> cachedFieldData;
  fill_cached_field_data(field, fieldIndex, cachedFieldData);
  return cachedFieldData;
}

void fill_cached_field_data(const std::vector<stk::mesh::FieldBase*>& fieldVec, std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData)
{
  size_t numFields = fieldVec.size();

  cachedFieldData.clear();

  for(unsigned i=0; i<numFields; i++) {
    const stk::mesh::FieldBase* field = fieldVec[i];
    cachedFieldData.push_back(get_cached_field_data(field, i));
  }
}

void fill_cached_const_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex, std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  if(nullptr == field) return;

  field_datatype_execute(*field,
    [&]<typename T>(const stk::mesh::FieldBase& /*fieldBase*/) {
      impl::internal_fill_cached_const_field_data<T>(field, fieldIndex, cachedFieldData);
    }
  );
}

void fill_cached_const_field_data(const stk::mesh::FieldBase* field, std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  fill_cached_const_field_data(field, 0u, cachedFieldData);
}

std::shared_ptr<CachedFieldDataBase> get_cached_const_field_data(const stk::mesh::FieldBase* field, const unsigned fieldIndex)
{
  std::shared_ptr<CachedFieldDataBase> cachedFieldData;
  fill_cached_const_field_data(field, fieldIndex, cachedFieldData);
  return cachedFieldData;
}

void fill_cached_const_field_data(const std::vector<stk::mesh::FieldBase*>& fieldVec, std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData)
{
  size_t numFields = fieldVec.size();

  cachedFieldData.clear();

  for(unsigned i=0; i<numFields; i++) {
    const stk::mesh::FieldBase* field = fieldVec[i];
    cachedFieldData.push_back(get_cached_const_field_data(field, i));
  }
}

void clear_cached_field_data(std::vector< std::shared_ptr<CachedFieldDataBase> > &cachedFieldData)
{
  cachedFieldData.clear();
}

void clear_cached_field_data(std::shared_ptr<CachedFieldDataBase>& cachedFieldData)
{
  cachedFieldData.reset();
}

} // namespace search
} // namespace stk



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

#ifndef TRANSFER_COPY_TRANSLATOR_HPP
#define TRANSFER_COPY_TRANSLATOR_HPP

#include <stk_util/stk_config.h>
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc
#include <functional>                       // for std::function

namespace stk {
namespace transfer {

struct DataTypeKey
{
  enum key_t : uint32_t {
      RANK_SHIFT = 28U
    , MAX_LENGTH = (1U << RANK_SHIFT) - 1U
    , TYPE_MASK = MAX_LENGTH
    , INVALID = ~0U
  };

  enum data_t {
    UNSIGNED_INTEGER,
    LONG_INTEGER,
    UNSIGNED_INTEGER_64,
    INTEGER,
    LONG_DOUBLE,
    DOUBLE,
    INVALID_TYPE
  };

  KOKKOS_FUNCTION
  static bool is_valid_data_type(data_t dataFieldType)
  {
    return dataFieldType >= UNSIGNED_INTEGER && dataFieldType <= DOUBLE;
  }

  KOKKOS_FUNCTION
  DataTypeKey()
    : m_value(INVALID)
  {}

  KOKKOS_FUNCTION
  DataTypeKey(uint32_t argKey)
    : m_value(static_cast<key_t>(argKey))
  {}

  KOKKOS_FUNCTION
  DataTypeKey(data_t argType, unsigned argLength)
    : m_value(static_cast<key_t>(static_cast<uint32_t>(argType) << RANK_SHIFT | argLength) )
  {
    STK_NGP_ThrowRequireMsg(argType <= static_cast<data_t>(255), "Error: given an out of range data type");
    STK_NGP_ThrowRequireMsg(argLength <= MAX_LENGTH, "Error: given an out of range length value");
  }

  KOKKOS_FUNCTION
  unsigned get_data_length() const { return m_value & TYPE_MASK; }

  KOKKOS_FUNCTION
  data_t get_data_type() const { return static_cast<data_t>(m_value >> RANK_SHIFT); }

  KOKKOS_FUNCTION
  bool is_valid() const { return m_value != INVALID; }

  KOKKOS_FUNCTION
  operator key_t() const { return m_value; }

  key_t m_value;
};

template<typename SRCTYPE, typename DESTTYPE>
inline void assign_src_type_to_dest_type(const void* srcAddr, void* destAddr, unsigned index) {
  ((DESTTYPE*)destAddr)[index] = ((SRCTYPE*)srcAddr)[index];
}

class TranslatorBase
{
 public:
  TranslatorBase() {}

  virtual void translate(const void* srcAddr, unsigned srcDataByteSize, DataTypeKey::data_t destType, void* destAddr, unsigned destDataByteSize) const = 0;
};

struct TranslatorInfo
{
  typedef std::function<void(const void*, void*, unsigned)> TranslateFunction;
  TranslateFunction translateFunction;
  size_t typeSize;
};

template<typename TYPE>
class DataTypeTranslator : public TranslatorBase
{

 public:
  DataTypeTranslator()
  {
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, unsigned>, sizeof(unsigned)});
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, int64_t>, sizeof(int64_t)});
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, uint64_t>, sizeof(uint64_t)});
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, int>, sizeof(int)});
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, long double>, sizeof(long double)});
    m_translators.push_back({assign_src_type_to_dest_type<TYPE, double>, sizeof(double)});
  }
  
  void translate(const void* srcAddr, unsigned srcDataByteSize, DataTypeKey::data_t destType, void* destAddr, unsigned destDataByteSize) const {
    STK_ThrowRequire(DataTypeKey::is_valid_data_type(destType));
    STK_ThrowRequire(srcDataByteSize % sizeof(TYPE) == 0);
    unsigned srcCount = srcDataByteSize / sizeof(TYPE);
    unsigned destCount = destDataByteSize / m_translators[destType].typeSize;

    unsigned count = std::min(srcCount, destCount);

    for(unsigned i = 0; i < count; i++) {
      m_translators[destType].translateFunction(srcAddr, destAddr, i);
    }
  }

 private:
  std::vector<TranslatorInfo> m_translators;
};

} }

#endif

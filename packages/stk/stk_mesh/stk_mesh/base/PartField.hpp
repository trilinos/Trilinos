/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_PartField_hpp
#define stk_mesh_PartField_hpp

#include <stk_mesh/base/Types.hpp>
#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class MetaData;

class PartFieldBase {
public:
  PartFieldBase(const MetaData* meta_ptr, unsigned partFieldIndex, unsigned itemsPerPart, unsigned bytesPerPart);

  virtual ~PartFieldBase() {
    for(size_t i=0; i<m_data_ptrs.size(); ++i) {
      delete [] m_data_ptrs[i];
    }
  }

  const MetaData* mesh_meta_data() const { return m_meta_data; }

  unsigned part_field_index() const { return m_index; }

  unsigned items_per_part() const { return m_items_per_part; }

  unsigned bytes_per_part() const { return m_bytes_per_part; }

  std::vector<char*>& char_data() { return m_data_ptrs; }
  const std::vector<char*>& char_data() const { return m_data_ptrs; }

private:
  const MetaData* m_meta_data;
  unsigned m_index;
  unsigned m_items_per_part;
  unsigned m_bytes_per_part;
protected:
  std::vector<char*> m_data_ptrs;
};


template< typename DataType >
class PartField : public PartFieldBase {
public:
  typedef DataType PartFieldDataType;

  PartField(const MetaData* meta_ptr, unsigned partFieldIndex, unsigned itemsPerPart)
  : PartFieldBase(meta_ptr, partFieldIndex, itemsPerPart, itemsPerPart*sizeof(DataType))
  {
    size_t bytes_per_item = bytes_per_part()/items_per_part();

    for(size_t i=0; i<m_data_ptrs.size(); ++i) {
      for(size_t j=0; j<items_per_part(); ++j) {
        char* ptr = m_data_ptrs[i] + j*bytes_per_item;
        new(ptr)PartFieldDataType;
      }
    }
  }

  virtual ~PartField() {
    size_t bytes_per_item = bytes_per_part()/items_per_part();

    for(size_t i=0; i<m_data_ptrs.size(); ++i) {
      for(size_t j=0; j<items_per_part(); ++j) {
        char* ptr = m_data_ptrs[i] + j*bytes_per_item;
        PartFieldDataType* item = reinterpret_cast<PartFieldDataType*>(ptr);
        item->~PartFieldDataType();
      }
    }
  }

  PartFieldDataType* data(unsigned part_ordinal)
  {
    char* char_ptr = m_data_ptrs[part_ordinal];
    return reinterpret_cast<PartFieldDataType*>(char_ptr);
  }
};

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_PartField_hpp */


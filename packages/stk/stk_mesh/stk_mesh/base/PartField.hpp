/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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


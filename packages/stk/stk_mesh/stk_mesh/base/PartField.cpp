
#include <stk_mesh/base/PartField.hpp>
#include <stk_mesh/base/MetaData.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

PartFieldBase::PartFieldBase(const MetaData* meta_ptr, unsigned partFieldIndex, unsigned itemsPerPart, unsigned bytesPerPart)
  : m_meta_data(meta_ptr), m_index(partFieldIndex), m_items_per_part(itemsPerPart), m_bytes_per_part(bytesPerPart)
{
  size_t num_parts = m_meta_data->get_parts().size();
  std::vector<char*>& char_ptr_vector = char_data();
  char_ptr_vector.resize(num_parts);
  for(size_t i=0; i<num_parts; ++i) {
    char_ptr_vector[i] = new char[bytes_per_part()];
  }
}

} // namespace mesh
} // namespace stk


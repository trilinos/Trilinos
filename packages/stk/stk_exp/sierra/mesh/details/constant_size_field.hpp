#ifndef SIERRA_SIERRA_MESH_FIELD_CONSTANT_SIZE_HPP
#define SIERRA_SIERRA_MESH_FIELD_CONSTANT_SIZE_HPP

#include <sierra/mesh/details/selector.hpp>
#include <sierra/mesh/details/bucket_location.hpp>
#include <sierra/mesh/field_traits.hpp>

#include <boost/foreach.hpp>

#include <vector>

namespace sierra {
namespace mesh {
namespace details {

/**
  constant_size_field provides a vector-valued field which has the
  same number of scalars at every entity on which it is defined.
*/
template<typename T, int length_per_entity=1>
class constant_size_field {
 public:
  enum { field_length_per_entity = length_per_entity };

  typedef T value_type;
  typedef T* pointer;

  constant_size_field()
   : m_bucket_offsets(),
     m_bucket_data_ptrs()
   , m_field_data()
  {}

  template <class Mesh>
  void update_from_mesh(const selector& select,
                        const Mesh & mesh)
  {
    typename Mesh::bucket_range bdrange = mesh.get_buckets();
    m_bucket_data_ptrs.assign(std::distance(bdrange.first, bdrange.second), static_cast<pointer>(NULL));
    m_bucket_offsets.assign(std::distance(bdrange.first, bdrange.second), 0);
    size_t num_field_data_scalars = 0;
    BOOST_FOREACH(bucket_key bd, bdrange) {
      if ( contains(select, bd, mesh) ) {
        num_field_data_scalars += mesh.num_entities(bd)*field_length_per_entity;
      }
    }

    m_field_data.assign(num_field_data_scalars, 0);

    size_t offset = 0;
    BOOST_FOREACH(bucket_key bd, bdrange) {
      if ( contains(select, bd, mesh) ) {
        m_bucket_data_ptrs[bd] = &m_field_data[offset];
        m_bucket_offsets[bd] = offset;
        offset += mesh.num_entities(bd)*field_length_per_entity;
      }
    }
  }

  int offset(bucket_location bl) const
  {
    return m_bucket_offsets[bl.bucket()] + bl.ordinal()*field_length_per_entity;
  }

  pointer operator[](bucket_location bl)
  {
    return m_bucket_data_ptrs[bl.bucket()]+bl.ordinal()*field_length_per_entity;
  }

  const pointer operator[](bucket_location bl) const
  {
    return m_bucket_data_ptrs[bl.bucket()]+bl.ordinal()*field_length_per_entity;
  }

  pointer operator[](bucket_key bd)
  {
    return m_bucket_data_ptrs[size_t(bd)];
  }

  const pointer operator[](bucket_key bd) const
  {
    return m_bucket_data_ptrs[size_t(bd)];
  }

//  std::pair<pointer,pointer> operator[](bucket_key bd)
//  {
//    pointer first = m_bucket_data_ptrs[size_t(bd)];
//    pointer last = m_bucket_data_ptrs[size_t(bd)+1];
//    return std::make_pair( first, last);
//  }
//
//  std::pair<const pointer, const pointer> operator[](bucket_key bd) const
//  {
//    const pointer first = m_bucket_data_ptrs[size_t(bd)];
//    const pointer last = m_bucket_data_ptrs[size_t(bd)+1];
//    return std::make_pair( first, last);
//  }

  std::vector<value_type>& field_data_flat_array()
  {
    return m_field_data;
  }

  const std::vector<value_type>& field_data_flat_array() const
  {
    return m_field_data;
  }

 private:
  std::vector<int> m_bucket_offsets;
  std::vector<pointer> m_bucket_data_ptrs;
  std::vector<value_type> m_field_data;
};//class constant_size_field

} // namespace details

//specialization of sierra::mesh::field_traits for FieldConstantSize:
template<>
struct field_traits<details::constant_size_field<double> > {
  typedef details::constant_size_field<double>::value_type data_type;
};

template<>
struct field_traits<details::constant_size_field<int> > {
  typedef details::constant_size_field<int>::value_type data_type;
};

} // namespace mesh
} // namespace sierra

#endif //SIERRA_SIERRA_MESH_DETAILS_FIELDCONSTANTSIZE_HPP

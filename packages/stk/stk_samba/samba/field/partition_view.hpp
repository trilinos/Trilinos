#ifndef SAMBA_SAMBA_MESH_FIELD_PARTITION_VIEW_HPP
#define SAMBA_SAMBA_MESH_FIELD_PARTITION_VIEW_HPP

namespace samba { namespace detail {

template < typename DataType
          ,unsigned Rank = 1u
          ,typename Enable = void
         >
struct field_partition_view;


template < typename DataType >
struct field_partition_view< DataType, 1>
{
  typedef DataType          data_type;
  typedef size_t            size_type;

  typedef data_type         value_type;
  typedef const data_type   const_value_type;

  typedef data_type &       reference;
  typedef data_type const & const_reference;


  const_reference operator[](partition_index i) const
  { return m_data[m_dimension*i.offset()()]; }

  const_reference operator[](partition_offset i) const
  { return m_data[m_dimension*i()]; }

  const_reference operator[](size_type i) const
  { return m_data[m_dimension*i]; }


  reference operator[](partition_index i)
  { return m_data[m_dimension*i.offset()()]; }

  reference operator[](partition_offset i)
  { return m_data[m_dimension*i()]; }

  reference operator[](size_type i)
  { return m_data[m_dimension*i]; }

  size_type dimension() const
  { return m_dimension; }

  data_type * m_data;
  size_type m_dimension;

};


//non-const field_partition_view
template < typename DataType >
struct field_partition_view< DataType, 2 >
{
  typedef DataType          data_type;
  typedef size_t            size_type;

  typedef data_type *       value_type;
  typedef data_type const * const_value_type;

  typedef data_type *       reference;
  typedef data_type const * const_reference;

  const_reference operator[](partition_index i) const
  { return m_data + m_dimension*i.offset()(); }

  const_reference operator[](partition_offset i) const
  { return m_data + m_dimension*i(); }

  const_reference operator[](size_type i) const
  { return m_data + m_dimension*i; }


  reference operator[](partition_index i)
  { return m_data + m_dimension*i.offset()(); }

  reference operator[](partition_offset i)
  { return m_data + m_dimension*i(); }

  reference operator[](size_type i)
  { return m_data + m_dimension*i; }

  size_type dimension() const
  { return m_dimension; }

  data_type * m_data;
  size_type m_dimension;

};

//TODO need specialization for row_storage


}} //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_FIELD_PARTITION_VIEW_HPP


#ifndef SAMBA_SAMBA_ENTITY_ENTITY_PROXY_HPP
#define SAMBA_SAMBA_ENTITY_ENTITY_PROXY_HPP

#include <samba/types.hpp>
#include <samba/mesh/convertable_to_partition_proxy.hpp>
#include <samba/mesh/convertable_to_entity_proxy.hpp>

namespace samba {

namespace detail {

class partition_impl;

} //namespace detail

/**
 * A handle to the the entity. Can become invalidated; further, this type
 * should never be stored; use it like an iterator.
 *
 * This allows for convenient and higher-performing queries on entities.
 */
class entity_proxy
{
 public:
  //*************************************************************************
  //connectivity accessors
  //*************************************************************************

  partition_index_range nodes() const;
  partition_index_range edges() const;
  partition_index_range faces() const;
  partition_index_range elements() const;

  partition_index_iterator begin_nodes() const;
  partition_index_iterator begin_edges() const;
  partition_index_iterator begin_faces() const;
  partition_index_iterator begin_elements() const;

  partition_index_iterator end_nodes() const;
  partition_index_iterator end_edges() const;
  partition_index_iterator end_faces() const;
  partition_index_iterator end_elements() const;

  size_t num_connectivity(entity_rank rank) const;
  size_t num_nodes() const;
  size_t num_edges() const;
  size_t num_faces() const;
  size_t num_elements() const;

  template <typename T>
  std::pair<const T*, const T*> connectivity(entity_rank rank) const;
  template <typename T>
  std::pair<const T*, const T*> nodes() const;
  template <typename T>
  std::pair<const T*, const T*> edges() const;
  template <typename T>
  std::pair<const T*, const T*> faces() const;
  template <typename T>
  std::pair<const T*, const T*> elements() const;

  template <typename T>
  const T* begin_connectivity(entity_rank rank) const;
  template <typename T>
  const T* begin_nodes() const;
  template <typename T>
  const T* begin_edges() const;
  template <typename T>
  const T* begin_faces() const;
  template <typename T>
  const T* begin_elements() const;

  template <typename T>
  const T* end_connectivity(entity_rank rank) const;
  template <typename T>
  const T* end_nodes() const;
  template <typename T>
  const T* end_edges() const;
  template <typename T>
  const T* end_faces() const;
  template <typename T>
  const T* end_elements() const;

  //*************************************************************************
  //entity querys
  //*************************************************************************

  detail::convertable_to_partition_proxy  partition() const
  { return m_partition; }

  entity_part_range parts() const;
  entity_key        key() const;
  partition_index descriptor() const;

  entity_topology topology() const
  { return key().topology(); }

  process_id process() const
  { return key().process(); }

  entity_local_id local_id() const
  { return key().local_id(); }

  entity_rank rank() const
  { return descriptor().rank(); }

  samba::partition_id partition_id() const
  { return descriptor().partition(); }

  partition_offset offset() const
  { return m_offset; }

  template <typename Index>
  Index get() const;

  //*************************************************************************
  //constructors
  //*************************************************************************
  entity_proxy()
    : m_partition(NULL)
    , m_offset(partition_offset::invalid())
  {}

  entity_proxy( detail::partition_impl const* arg_partition
                ,partition_offset arg_offset
                )
    : m_partition(arg_partition)
    , m_offset(arg_offset)
  {
    BOOST_ASSERT_MSG(m_offset() < m_partition->size(),
                     (debug_message() << "partition_offset " << m_offset() << " out-of-bounds, max is: " << m_partition->size() - 1));
  }

  entity_proxy(detail::convertable_to_entity_proxy const& arg_proxy)
    : m_partition(arg_proxy.m_partition)
    , m_offset(arg_proxy.m_offset)
  {
    BOOST_ASSERT_MSG(m_offset() < m_partition->size(),
                     (debug_message() << "partition_offset " << m_offset() << " out-of-bounds, max is: " << m_partition->size() - 1));
  }

 private:
  detail::partition_impl const* m_partition;
  partition_offset                 m_offset;
};

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_ENTITY_IMPL_HPP


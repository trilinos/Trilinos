#ifndef SAMBA_SAMBA_MESH_PARTITION_PROXY_HPP
#define SAMBA_SAMBA_MESH_PARTITION_PROXY_HPP

#include <samba/set_expression.hpp>

namespace samba {

/**
 * Can become invalidated; further, this type should never be stored;
 * use it like an iterator. Simply provides a convenient handle to a
 * partition.
 */
class partition_proxy
{
 public:
  typedef detail::partition_impl::iterator iterator;
  typedef iterator const_iterator;

  //*************************************************************************
  //constructors
  //*************************************************************************
  partition_proxy()
    : m_partition(NULL)
  {}

  partition_proxy( detail::partition_impl const* arg_partition )
    : m_partition(arg_partition)
  {}

  partition_proxy( detail::convertable_to_partition_proxy arg_partition )
    : m_partition(arg_partition.m_partition)
  {}

  //*************************************************************************
  //connectivity accessors
  //*************************************************************************

  partition_index_range nodes(partition_offset offset) const;
  partition_index_range edges(partition_offset offset) const;
  partition_index_range faces(partition_offset offset) const;
  partition_index_range elements(partition_offset offset) const;

  partition_index_iterator begin_nodes(partition_offset offset) const;
  partition_index_iterator begin_edges(partition_offset offset) const;
  partition_index_iterator begin_faces(partition_offset offset) const;
  partition_index_iterator begin_elements(partition_offset offset) const;

  partition_index_iterator end_nodes(partition_offset offset) const;
  partition_index_iterator end_edges(partition_offset offset) const;
  partition_index_iterator end_faces(partition_offset offset) const;
  partition_index_iterator end_elements(partition_offset offset) const;

  // API below is always valid

  size_t num_connectivity(entity_rank rank, partition_offset offset) const;
  size_t num_nodes(partition_offset offset) const;
  size_t num_edges(partition_offset offset) const;
  size_t num_faces(partition_offset offset) const;
  size_t num_elements(partition_offset offset) const;

  template <typename T>
  std::pair<const T*, const T*> connectivity(entity_rank rank, partition_offset offset) const;
  template <typename T>
  std::pair<const T*, const T*> nodes(partition_offset offset) const;
  template <typename T>
  std::pair<const T*, const T*> edges(partition_offset offset) const;
  template <typename T>
  std::pair<const T*, const T*> faces(partition_offset offset) const;
  template <typename T>
  std::pair<const T*, const T*> elements(partition_offset offset) const;

  template <typename T>
  const T* begin_connectivity(entity_rank rank, partition_offset offset) const;
  template <typename T>
  const T* begin_nodes(partition_offset offset) const;
  template <typename T>
  const T* begin_edges(partition_offset offset) const;
  template <typename T>
  const T* begin_faces(partition_offset offset) const;
  template <typename T>
  const T* begin_elements(partition_offset offset) const;

  template <typename T>
  const T* end_connectivity(entity_rank rank, partition_offset offset) const;
  template <typename T>
  const T* end_nodes(partition_offset offset) const;
  template <typename T>
  const T* end_edges(partition_offset offset) const;
  template <typename T>
  const T* end_faces(partition_offset offset) const;
  template <typename T>
  const T* end_elements(partition_offset offset) const;

  //*************************************************************************
  //partition querys
  //*************************************************************************

  samba::spatial_dimension spatial_dimension() const;
  entity_rank rank() const;
  samba::partition_id partition_id() const;
  entity_topology topology() const;
  entity_part_range parts() const;

  size_t size() const;
  bool empty() const;

  entity_proxy first() const;
  entity_proxy last() const;

  const_iterator begin() const;
  const_iterator end() const;

  samba::connectivity_kind connectivity_kind(entity_rank r) const;

  //*************************************************************************
  //entity querys
  //*************************************************************************

  entity_proxy operator[](partition_index d) const;
  entity_proxy operator[](partition_offset o) const;
  entity_proxy operator[](unsigned i) const;

  entity_key        key(partition_offset offset) const;
  partition_index descriptor(partition_offset offset) const;

  template <typename RankIndex>
  RankIndex rank_index(partition_offset offset) const;

 private:
  detail::partition_impl const*  m_partition;
};

inline bool contains(set_expression const& set, partition_proxy partition )
{ return contains(set, partition.parts()); }

inline bool contains(set_expression const& set, detail::convertable_to_partition_proxy partition )
{ return contains(set, partition_proxy(partition).parts()); }

} //namespace samba

#endif //SAMBA_SAMBA_MESH_ENTITY_ENTITY_SUBSET_HPP

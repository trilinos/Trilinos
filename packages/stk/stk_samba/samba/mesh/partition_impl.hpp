#ifndef SAMBA_SAMBA_MESH_PARTITION_IMPL_HPP
#define SAMBA_SAMBA_MESH_PARTITION_IMPL_HPP

#include <samba/types.hpp>
#include <samba/mesh/convertable_to_partition_proxy.hpp>
#include <samba/mesh/convertable_to_entity_proxy.hpp>
#include <samba/mesh/partition_connectivity.hpp>

namespace samba { namespace detail {

class partition_iterator;
class mesh_impl;

/**
 * partition_proxy represents a collection of entities that are homogeneous
 * (have the same parts, entity_rank, and topology). This class manages
 * the descriptors, keys, and connectivity for the entities it contains.
 *
 * This is similar to the "bucket" concept in our older meshes, except no
 * field data is being managed here.
 */
class partition_impl
{
 public:
  typedef partition_iterator iterator;
  typedef iterator           const_iterator;

  //*************************************************************************
  //constructors
  //*************************************************************************

  partition_impl( mesh_impl * arg_mesh
                  ,partition_id arg_partition
                  ,entity_part_vector const& arg_parts
                  );

  // partition_index API not valid during modification

  // partition_index_* changes to node_index_* (and same for edges, faces, elements)
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

  // API below is always valid. Clients should use the iterator typedefs in types.hpp
  // to hold the return values.

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
  partition_id partition() const;
  entity_topology topology() const;
  entity_rank rank() const;
  entity_part_vector const& parts() const;

  size_t size() const;
  bool empty() const;

  convertable_to_entity_proxy first() const;
  convertable_to_entity_proxy last() const;

  const_iterator begin() const;
  const_iterator end() const;

  samba::connectivity_kind connectivity_kind(entity_rank r) const;

  template <typename EntityRank>
  samba::connectivity_kind connectivity_kind() const;

  uint32_t starting_rank_count() const
  { return m_starting_rank_count; }

  //*************************************************************************
  //entity querys
  //*************************************************************************

  convertable_to_entity_proxy operator[](partition_index d) const;
  convertable_to_entity_proxy operator[](partition_offset o) const;
  convertable_to_entity_proxy operator[](unsigned i) const;

  entity_key        key(partition_offset offset) const;
  partition_index   descriptor(partition_offset offset) const;

  template <typename RankIndex>
  RankIndex rank_index(partition_offset offset) const;

  //*************************************************************************
  //modification api
  //*************************************************************************

  void update_partition_id(partition_id descriptor)
  { m_partition = descriptor; }

  // When you add
  //   case 1: connectivity_kind == invalid -> ASSERT_MSG
  //   case 2: connectivity_kind == fixed -> Sets connectivity (memory should already be allocated)
  //   case 3: connectivity_kind == dynamic -> Add column to row_storage
  void add_connectivity   (partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation=connectivity_orientation::invalid());

  // When you remove: same cases as above, just inverse
  void remove_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal);

  void swap(partition_offset first, partition_offset second);

  // Call when modification cycle is completed. Inflates partition_index vectors
  void end_modification_impl();

  void begin_modification_impl();

  void add_entities_impl(size_t how_many);

  void remove_entities_impl(size_t how_many);

  void move_entities_impl(partition_impl& to, size_t how_many);

  void set_starting_rank_count(uint32_t count)
  { m_starting_rank_count = count; }

 private:

  void check_invariant_helper() const;

  // Data Members:
  mesh_impl * m_mesh;
  partition_id m_partition;
  entity_part_vector m_parts;

  size_t m_num_entities;
  uint32_t m_starting_rank_count;
  bool m_in_modification_cycle;

  samba::connectivity_kind m_node_kind;
  samba::connectivity_kind m_edge_kind;
  samba::connectivity_kind m_face_kind;
  samba::connectivity_kind m_element_kind;

  partition_connectivity<entity_rank::node_type, samba::connectivity_kind::fixed_type> m_fixed_node_connectivity; // fixed connectivity to nodes
  partition_connectivity<entity_rank::edge_type, samba::connectivity_kind::fixed_type> m_fixed_edge_connectivity; // fixed connectivity to edges
  partition_connectivity<entity_rank::face_type, samba::connectivity_kind::fixed_type> m_fixed_face_connectivity; // fixed connectivity to faces
  partition_connectivity<entity_rank::element_type, samba::connectivity_kind::fixed_type> m_fixed_element_connectivity; // fixed connectivity to elements

  partition_connectivity<entity_rank::node_type, samba::connectivity_kind::dynamic_type> m_dynamic_node_connectivity; // dynamic connectivity to nodes
  partition_connectivity<entity_rank::edge_type, samba::connectivity_kind::dynamic_type> m_dynamic_edge_connectivity; // dynamic connectivity to edges
  partition_connectivity<entity_rank::face_type, samba::connectivity_kind::dynamic_type> m_dynamic_face_connectivity; // dynamic connectivity to faces
  partition_connectivity<entity_rank::element_type, samba::connectivity_kind::dynamic_type> m_dynamic_element_connectivity; // dynamic connectivity to elements
};

}} //namespace samba::detail

#endif // SAMBA_SAMBA_MESH_PARTITION_IMPL_HPP

#ifndef SAMBA_SAMBA_MESH_PARTITION_CONNECTIVITY_HPP
#define SAMBA_SAMBA_MESH_PARTITION_CONNECTIVITY_HPP

#include <samba/types.hpp>
#include <boost/assert.hpp>

namespace samba { namespace detail {

template< typename ToEntityRank, typename ConnectivityKind >
class partition_connectivity;
//
//{
//public:
//
//  typedef partition_connectivity<ToEntityRank,ConnectivityKind> self_type;
//
//  partition_connectivity(entity_topology from_topology, spatial_dimension arg_spatial_dimension);
//
//  // Query API
//  // Vaild Index:
//  //   entity_key
//  //   partition_index
//  //   connectivity_ordinal
//  //   connectivity_orientation
//  //   rank_index
//  template< typename Index >
//  std::pair<const Index*, const Index *> range(partition_offset offset) const;
//
//  template< typename Index >
//  const Index * begin(partition_offset offset) const;
//
//  template< typename Index >
//  const Index * end(partition_offset offset) const;
//
//  size_t num_connectivity(partition_offset offset) const;
//
//  // Modification API
//  void add_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation);
//
//  void remove_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal);
//
//  void swap(partition_offset first, partition_offset second);
//
//  void begin_modification();
//
//  void end_modification(mesh_impl * mesh);
//
//  void add_entities(size_t how_many);
//
//  void remove_entities(size_t how_many);
//
//  void move_entities(self_type & to, size_t how_many);
//
//};

}} //namespace samba::detail

#endif // SAMBA_SAMBA_MESH_PARTITION_CONNECTIVITY_HPP

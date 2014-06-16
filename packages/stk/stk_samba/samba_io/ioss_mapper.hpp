#ifndef SAMBA_SAMBA_IO_IOSS_MAPPER_HPP
#define SAMBA_SAMBA_IO_IOSS_MAPPER_HPP

#include <Ioss_SubSystem.h>

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba/rank_field.hpp>
#include <samba_io/detail/sidesets_info.hpp>

namespace samba {
namespace io {

typedef mesh::entity_key_vector entity_key_vector;

typedef rank_field<entity_rank::node_type, double, spatial_dimension_functor> coordinates_field_type;

typedef field<int64_t, scalar_functor>            entity_key_to_ioss_id_field;

typedef boost::unordered_map<entity_block_key, Ioss::EntityBlock *>  entity_block_to_iossblock_map;

typedef boost::unordered_map<entity_block_key, Ioss::EntityBlock *>  entity_block_to_iossblock_map;


/// Data structure for data that augments a samba::mesh to facilitate writing out a
/// model (e.g., Exodus) file.
struct ioss_mapper
{
  typedef boost::shared_ptr<ioss_mapper> Ptr;

  ioss_mapper (mesh mesh_arg)
    : m_to_ioss_idx_in_block(mesh_arg, entity_state::universe(), -1)
    , m_to_ioss_local_id(mesh_arg, entity_state::universe(), -1)
    , m_to_ioss_global_id(mesh_arg, entity_state::universe(), -1)
    , m_node_coordinates(mesh_arg)
  {
    // Nada.
  }

  /// Map entity keys to index in the block in the Exodus file that created the entity.
  entity_key_to_ioss_id_field                                          m_to_ioss_idx_in_block;

  /// Map entity keys to local IDs in an Exodus file.  We can use this to preserve
  /// the Exodus local IDs of nodes as long as no nodes are destroyed.  Not sure
  /// whether to keep this.
  entity_key_to_ioss_id_field                                          m_to_ioss_local_id;

  /// Map entity keys to global Exodus IDs.
  entity_key_to_ioss_id_field                                          m_to_ioss_global_id;

  /// Node coordinates.
  coordinates_field_type                                               m_node_coordinates;

  /// Samba entity blocks (sets of entities) that map to "node blocks" in Ioss; should only be
  /// one for an Exodus file.
  std::vector<entity_block_key>                                        m_node_blocks;

  /// Needed to enable entities added to a block after it is read in to be written out.
  std::map<entity_block_key, int>                                      m_new_idx_in_block;

  /// Structure of sidesets and side blocks.
  std::vector<SidesetPtr>                                              m_sidesets_info;
};


} // namespace samba::io

} // namespace samba

#endif

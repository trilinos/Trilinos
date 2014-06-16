#ifndef SAMBA_SAMBA_MESH_MESH_SIGNALS_HPP
#define SAMBA_SAMBA_MESH_MESH_SIGNALS_HPP

#include <samba/types.hpp>

#include <boost/signals2.hpp>

namespace samba {

//*************************************************************************
//signals - A signal/slot (similar to callback concept) system that allows
//          implementors of fields to be informed about changes in mesh
//          structure.
//
//          e.g.
//          max_local_id_signal.connect(boost::bind(&<func>, this, _1, _2));
//          ...
//          void <func>(entity_rank arg_rank, size_t max_local_id) {...}
//*************************************************************************
struct mesh_signals
{
  /**
   * signal when user creates an entity_block
   */
  typedef boost::signals2::signal<void ( entity_block_key /*set*/ ) > add_entity_block_signal_type;

  add_entity_block_signal_type add_entity_block_signal;


  /**
   * signal when entity keys are added to the mesh
   */
  typedef boost::signals2::signal<void ( entity_key_interval )> add_entity_keys_signal_type;

  add_entity_keys_signal_type add_entity_keys_signal;

  /**
   * signal when entity keys are removed from the mesh
   */
  typedef boost::signals2::signal<void ( entity_key_iterator /*first*/
                              ,entity_key_iterator /*last*/)
                       > remove_entity_keys_signal_type;

  remove_entity_keys_signal_type remove_entity_keys_signal;

  /**
   * signal when a partition is added to the mesh
   */
  typedef boost::signals2::signal<void ( partition_id /*partition*/) > add_partition_signal_type;

  add_partition_signal_type add_partition_signal;

  /**
   * signal when partitions are reordered
   */
  typedef boost::signals2::signal<void ( std::vector<partition_id> const & /*order*/ ) > reorder_partitions_signal_type;

  reorder_partitions_signal_type reorder_partitions_signal;

  /**
   * signal when entities are added to a partition
   */
  typedef boost::signals2::signal<void ( partition_id /*partition*/
                              ,size_t /*how_many*/ )
                       > add_entities_signal_type;

  add_entities_signal_type add_entities_signal;

  /**
   * signal when entities are moved from one partition to another
   */
  typedef boost::signals2::signal<void ( partition_id /*partition_from*/
                              ,partition_id /*partition_to*/
                              ,size_t /*how_many*/ )
                       > move_entities_signal_type;

  move_entities_signal_type move_entities_signal;

  /**
   * signal when entities are remove from a partition
   */
  typedef boost::signals2::signal<void ( partition_id /*partition_from*/
                              ,size_t /*how_many*/ )
                       > remove_entities_signal_type;

  remove_entities_signal_type remove_entities_signal;

  /**
   * signal when entities are reordered within a partition
   */
  typedef boost::signals2::signal<void ( partition_id /*partition*/
                              ,std::vector<partition_offset> const & /*order*/)
                       > reorder_partition_signal_type;

  reorder_partition_signal_type reorder_partition_signal;

  /**
   * signal when entities are swapped within a partition
   */
  typedef boost::signals2::signal<void ( partition_id /*partition*/
                              ,partition_offset /*from*/
                              ,partition_offset /*to*/)
                       > swap_offset_signal_type;

  swap_offset_signal_type swap_offset_signal;

  /**
   * signal when modification of the mesh begins
   */
  typedef boost::signals2::signal<void ()> begin_modification_signal_type;

  begin_modification_signal_type begin_modification_signal;

  /**
   * signal when modification of the mesh ends
   */
  typedef boost::signals2::signal<void ()> end_modification_signal_type;

  end_modification_signal_type end_modification_signal;
};

} //namespace samba

#endif // SAMBA_SAMBA_MESH_MESH_SIGNALS_HPP

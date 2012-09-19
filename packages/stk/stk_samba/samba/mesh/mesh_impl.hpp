#ifndef SAMBA_SAMBA_MESH_MESH_IMPL_HPP
#define SAMBA_SAMBA_MESH_MESH_IMPL_HPP

#include <samba/rank_index.hpp>
#include <samba/types.hpp>
#include <samba/set_expression.hpp>

#include <samba/mesh/convertable_to_partition_proxy.hpp>
#include <samba/mesh/convertable_to_entity_proxy.hpp>
#include <samba/mesh/mesh_signals.hpp>
#include <samba/mesh/query_impl.hpp>
#include <samba/mesh/conversion_impl.hpp>

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <bitset>

namespace samba { namespace detail {

class partition_impl;

/**
 * Implementation of mesh interface.
 */
class mesh_impl : public query_impl<mesh_impl>, public conversion_impl<mesh_impl>
{
public:
  typedef boost::shared_ptr<partition_impl> partition_handle;

  //*************************************************************************
  //constructor
  //*************************************************************************

  mesh_impl(  samba::connectivity_map const& arg_connectivity_map = samba::connectivity_map::default_map()
             ,process_id arg_process = process_id::invalid()
           );

  //*************************************************************************
  //signals
  //*************************************************************************
  mesh_signals & signals();

  mesh_signals const& signals() const;

  //*************************************************************************
  //Queries
  //*************************************************************************


  //*************************************************************************
  //modification
  //*************************************************************************

  void add_set_expression_by_name(const std::string& name, set_expression expr);

  bool begin_modification();

  bool end_modification();

  entity_block_key add_entity_block(std::string const& name, entity_rank rank, bool inducable);

  template <typename EntitySetIterator>
  entity_key_interval add_entities( entity_topology topology
                                    ,size_t number
                                    ,EntitySetIterator first
                                    ,EntitySetIterator last
                                    ,bool is_owned = true
                                    );

  template <typename EntityKeyIterator>
  void remove_entities( EntityKeyIterator first
                        ,EntityKeyIterator last
                        );

  template <typename EntityKeyIterator, typename EntitySetIterator, typename RemoveEntitySetIterator>
  void move_entities( EntityKeyIterator first_key
                      ,EntityKeyIterator last_key
                      ,EntitySetIterator first_add_set
                      ,EntitySetIterator last_add_set
                      ,RemoveEntitySetIterator first_remove_set
                      ,RemoveEntitySetIterator last_remove_set
                      );


  void add_connectivity( entity_key from
                     ,entity_key to
                     ,connectivity_ordinal ordinal
                     ,connectivity_orientation orientation = connectivity_orientation::invalid()
                     );

  void remove_connectivity( entity_key from
                        ,entity_key to
                        ,connectivity_ordinal ordinal
                        );

  template <typename ForwardEntityKeyIterator>
  void add_non_owned_keys( std::vector<entity_key> & local_keys
                          ,ForwardEntityKeyIterator first
                          ,ForwardEntityKeyIterator last
                        );

private:

  static const size_t MAX_INDUCED_PARTS = 64;

  typedef boost::unordered_map<entity_part_vector,partition_id>  part_vector_to_partition_map;
  typedef std::map<std::string, set_expression>                  expr_db_type;
  typedef std::bitset<MAX_INDUCED_PARTS>                         induced_part_set;

  typedef boost::unordered_map<entity_key,entity_key>            distributed_key_map;

  void add_topology_helper(entity_topology topology);

  void max_local_id_changed_helper( entity_topology topology
                                    ,entity_local_id max_id
                                    );

  bool has_been_deleted_helper(entity_key key) const;

  entity_key_interval create_keys_helper(entity_topology topology, size_t how_many);

  interval<entity_local_id> get_local_ids_helper(entity_topology topology, size_t how_many);

  partition_id get_partition_non_const( entity_part_vector const& parts );

  partition_id create_partition( entity_part_vector const& parts );

  template <typename Iterator>
  void move_entities_helper( partition_id from
                             ,partition_id to
                             ,Iterator first
                             ,Iterator last
                             ,bool during_end_modification=false
                             );

  template <typename Iterator>
  void remove_entities_helper( partition_id partition
                               ,Iterator first
                               ,Iterator last
                               );

  void remove_relations_to_deleted_entities_helper();

  void apply_induced_part_membership_helper();

  void compress_data_structure_helper(std::vector<partition_id>& new_partition_order);

  void populate_rank_index_helper();

  void add_induced_parts_helper(entity_part_vector const& parent_parts,
                                const size_t first_induced_idx,
                                induced_part_set& induced_parts);

  size_t find_induced_parts_helper(entity_part_vector const& orig_parts) const;

  void induce_parts_through_downward_connectivity_helper(partition_handle                              p_handle,
                                                         std::vector< std::vector<induced_part_set> >& key_to_induced_part_map);

  template <typename RelationIterator>
  void induce_parts_helper(entity_part_vector const&                  parent_parts,
                           const size_t                               first_inducable_idx,
                           RelationIterator                           begin,
                           RelationIterator                           end,
                           std::vector< std::vector<induced_part_set> >& key_to_induced_part_map);

  void validate_helper(partition_id p) const;
  void validate_helper(partition_index d) const;
  void validate_helper(entity_key k) const;
  void validate_helper(entity_topology t) const;
  void validate_helper(entity_block_key k) const;
  void validate_helper(entity_rank r) const;
  void validate_helper(entity_part p) const;
  void validate_helper(entity_state s) const;
  void validate_helper(node_index i) const;
  void validate_helper(edge_index i) const;
  void validate_helper(face_index i) const;
  void validate_helper(element_index i) const;

  template <typename IteratorT>
  void validate_helper(IteratorT first, IteratorT last) const;

  template <typename EntityRank>
  const std::vector<partition_index>& get_rank_index_map() const;

  struct entity_block_spec
  {
    std::string name;
    entity_rank rank;
    bool inducable;
    int  inducable_index;
  };

  struct entity_block_match_name
  {
    entity_block_match_name(std::string const& name) : m_name(name) {}

    bool operator()(entity_block_spec const& spec)
    {
      return spec.name == m_name;
    }

    std::string const& m_name;
  };

  // MEMBERS

  process_id m_process;

  bool m_is_in_modification_cycle;

  samba::connectivity_map m_connectivity_map;

  // indexed by entity_block_key
  std::vector<entity_block_spec> m_set_specs;

  // indexed by inducable_index within entity_block_spec
  std::vector<entity_block_key> m_inducable_parts;

  // entity_part_vector -> partition_id
  part_vector_to_partition_map m_parts_to_partition;

  // indexed by partition_id
  std::vector<partition_handle> m_partitions;

  // indexed by entity_topology, returns interval set of local-ids
  std::vector<interval_set<entity_local_id> > m_local_ids;

  // represents deleted keys in the current modification cycle
  interval_set<entity_key> m_deleted_keys;

  // maps non-local keys to partition_indices
  boost::unordered_map<entity_key, partition_index> m_non_local_key_to_descriptor_map;

  // indexed by [entity_key.topology()][entity_key.local_id()], returns an partition_index, local keys only
  std::vector<std::vector<partition_index> > m_key_to_descriptor_map;

  // indexed by [partition_index.partition()][partition_index.offset()], returns an entity_key
  std::vector<std::vector<entity_key> > m_descriptor_to_key_map;

  mesh_signals m_signals;

  // set_expression database
  expr_db_type m_expr_db;

  // indexed by node_index
  std::vector<partition_index> m_node_index_to_partition_index_map;

  // indexed by edge_index
  std::vector<partition_index> m_edge_index_to_partition_index_map;

  // indexed by face_index
  std::vector<partition_index> m_face_index_to_partition_index_map;

  // indexed by element_index
  std::vector<partition_index> m_element_index_to_partition_index_map;

  friend class query_impl<mesh_impl>;
  friend class conversion_impl<mesh_impl>;
  friend class partition_impl;
};

}} //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_MESH_IMPL_HPP

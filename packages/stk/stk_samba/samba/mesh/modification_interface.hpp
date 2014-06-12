#ifndef SAMBA_SAMBA_MESH_MODIFICATION_INTERFACE_HPP
#define SAMBA_SAMBA_MESH_MODIFICATION_INTERFACE_HPP

namespace samba {

/**
 * The interface class for this package. Clients will interact with this class
 * to perform mesh operations.
 */
template <typename MeshHandle>
class modification_interface
{
 public:

  void add_set_expression_by_name(const std::string& name, set_expression const & expr);

  /**
   * Begin a modification cycle. Returns true if a modification cycle was started
   * (I.E. returns false if user was already in a modification cycle).
   *
   * You must be in a modification cycle in order to make changes to the mesh.
   * Warning: changing the mesh can invalidate descriptors, so during a modification
   * cycle, use keys instead of descriptors.
   */
  bool begin_modification();

  /**
   * End a modification cycle. Returns true if a modification cycle was ended
   * (I.E. returns false if there was no modification cycle to end).
   *
   * Ending a modification cycle will...
   *   - optimize mesh internal data structures,
   *   - induce part membership based on connectivity
   *
   * Descriptors will be stable outside of modification
   * cycles and are faster than keys so should be preferred.
   */
  bool end_modification();

  /**
   * Add a new user-controlled set of entities.
   *
   * Blocks can be associated with a rank; this restricts what entities can be
   * added to this block.
   *
   * Can be set to "inducable", meaning that, upon commit,
   * entities in this set will induce lower-ranked entities for which
   * it has a connectivity into this set. Induced sets are transitive.
   * Entities can be induced into a block with a rank even if they are
   * not of the same rank.
   */
  entity_block_key add_entity_block(std::string const& name, entity_rank rank=entity_rank::invalid(), bool inducable=false);

  template <typename EntitySetIterator>
  entity_key_interval add_entities( entity_topology topology
                                    ,size_t how_many
                                    ,EntitySetIterator first
                                    ,EntitySetIterator last
                                    );

  entity_key_interval add_entities( entity_topology topology
                                    ,size_t how_many
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

  template <typename EntityKeyIterator, typename EntitySetIterator>
  void move_entities( EntityKeyIterator first_key
                      ,EntityKeyIterator last_key
                      ,EntitySetIterator first_add_set
                      ,EntitySetIterator last_add_set
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

 protected:
  ~modification_interface() {}

 private:
  detail::mesh_impl& mesh_impl()
  { return *(static_cast<MeshHandle*>(this)->m_mesh_impl); }
};

} //namespace samba

#endif // SAMBA_SAMBA_MESH_MODIFICATION_INTERFACE_HPP

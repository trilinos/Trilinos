namespace samba {

template <typename MeshHandle>
inline
void modification_interface<MeshHandle>::add_set_expression_by_name(const std::string& name, set_expression const & expr)
{ return mesh_impl().add_set_expression_by_name(name, expr); }

/**
 * Begin a modification cycle. Returns true if a modification cycle was started
 * (I.E. returns false if user was already in a modification cycle).
 *
 * You must be in a modification cycle in order to make changes to the mesh.
 * Warning: changing the mesh can invalidate descriptors, so during a modification
 * cycle, use keys instead of descriptors.
 */
template <typename MeshHandle>
inline
bool modification_interface<MeshHandle>::begin_modification()
{ return mesh_impl().begin_modification(); }

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
template <typename MeshHandle>
inline
bool modification_interface<MeshHandle>::end_modification()
{ return mesh_impl().end_modification(); }

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
template <typename MeshHandle>
inline
entity_block_key modification_interface<MeshHandle>::add_entity_block(std::string const& name,
                                                                      entity_rank rank,
                                                                      bool inducable)
{ return mesh_impl().add_entity_block(name, rank, inducable); }


template <typename MeshHandle>
template <typename EntitySetIterator>
inline
entity_key_interval modification_interface<MeshHandle>::add_entities( entity_topology topology
                                                                      ,size_t how_many
                                                                      ,EntitySetIterator first
                                                                      ,EntitySetIterator last
                                                                      )
{
  return mesh_impl().add_entities(topology,how_many,first,last);
}

template <typename MeshHandle>
inline
entity_key_interval modification_interface<MeshHandle>::add_entities( entity_topology topology
                                                                      ,size_t how_many
                                                                      )
{
  entity_block_key none = entity_block_key::invalid();
  return add_entities(topology,how_many,&none,&none);
}

template <typename MeshHandle>
template <typename EntityKeyIterator>
inline
void modification_interface<MeshHandle>::remove_entities( EntityKeyIterator first
                                                          ,EntityKeyIterator last
                                                          )
{
  mesh_impl().remove_entities(first,last);
}

template <typename MeshHandle>
template <typename EntityKeyIterator, typename EntitySetIterator, typename RemoveEntitySetIterator>
inline
void modification_interface<MeshHandle>::move_entities( EntityKeyIterator first_key
                                                        ,EntityKeyIterator last_key
                                                        ,EntitySetIterator first_add_set
                                                        ,EntitySetIterator last_add_set
                                                        ,RemoveEntitySetIterator first_remove_set
                                                        ,RemoveEntitySetIterator last_remove_set
                                                        )
{
  mesh_impl().move_entities( first_key,last_key
                             ,first_add_set, last_add_set
                             ,first_remove_set, last_remove_set);
}

template <typename MeshHandle>
template <typename EntityKeyIterator, typename EntitySetIterator>
inline
void modification_interface<MeshHandle>::move_entities( EntityKeyIterator first_key
                                                        ,EntityKeyIterator last_key
                                                        ,EntitySetIterator first_add_set
                                                        ,EntitySetIterator last_add_set
                                                        )
{
  entity_block_key none;
  mesh_impl().move_entities( first_key,last_key
                             ,first_add_set, last_add_set
                             ,&none, &none);
}

template <typename MeshHandle>
inline
void modification_interface<MeshHandle>::add_connectivity( entity_key from
                                                           ,entity_key to
                                                           ,connectivity_ordinal ordinal
                                                           ,connectivity_orientation orientation
                                                           )
{
  return mesh_impl().add_connectivity(from,to,ordinal,orientation);
}

template <typename MeshHandle>
inline
void modification_interface<MeshHandle>::remove_connectivity( entity_key from
                                                              ,entity_key to
                                                              ,connectivity_ordinal ordinal
                                                              )
{
  return mesh_impl().remove_connectivity(from,to,ordinal);
}

} // namespace samba

/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

namespace stk {
namespace mesh {

namespace impl {

inline
unsigned get_ordinal(const Part* part)
{ return part->mesh_meta_data_ordinal(); }

inline
unsigned get_ordinal(unsigned ord)
{ return ord; }

inline
const Part& get_part(const Part* part, MetaData& meta)
{ return *part; }

inline
const Part& get_part(unsigned ord, MetaData& meta)
{ return *meta.get_parts()[ord]; }

template <typename PartT>
inline
void filter_out( std::vector<unsigned> & vec ,
                 const std::vector<PartT> & parts ,
                 std::vector<PartT> & removed )
{
  std::vector<unsigned>::iterator i , j ;
  i = j = vec.begin();

  typename std::vector<PartT>::const_iterator ip = parts.begin() ;

  while ( j != vec.end() && ip != parts.end() ) {
    if      ( get_ordinal(*ip) < *j ) { ++ip ; }
    else if ( *j < get_ordinal(*ip) ) { *i = *j ; ++i ; ++j ; }
    else {
      removed.push_back( *ip );
      ++j ;
      ++ip ;
    }
  }

  if ( i != j ) { vec.erase( i , j ); }
}

template <typename PartT>
inline
void merge_in( std::vector<unsigned> & vec , const std::vector<PartT> & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  typename std::vector<PartT>::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = get_ordinal(*ip);

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      // Now have: ord == *i
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = get_ordinal(*ip);
    vec.push_back( ord );
  }
}

} // namespace impl

template<typename AddIterator, typename RemoveIterator>
inline
void BulkData::change_entity_parts( Entity entity,
                                    AddIterator begin_add_parts, AddIterator end_add_parts,
                                    RemoveIterator begin_remove_parts, RemoveIterator end_remove_parts,
                                    bool always_propagate_internal_changes)
{
  TraceIfWatching("stk::mesh::BulkData::change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity);

  require_ok_to_modify();

  require_entity_owner( entity , m_parallel_rank );

  const EntityRank ent_rank = entity_rank(entity);
  const EntityRank undef_rank  = InvalidEntityRank;

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  // most parts will at least have universal and topology part as supersets
  const unsigned expected_min_num_supersets = 2;

  OrdinalVector a_parts;
  a_parts.reserve( std::distance(begin_add_parts, end_add_parts) * (expected_min_num_supersets + 1) );
  for(AddIterator add_iter=begin_add_parts; add_iter!=end_add_parts; ++add_iter) {
#ifdef FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    ThrowErrorMsgIf(entity_rank == stk::topology::ELEMENT_RANK && **add_iter == mesh_meta_data().globally_shared_part(), "FMWK_NO_GLOBALLY_SHARED_ELEMENTS  Error in BulkData::change_entity_parts, trying to make an element globally shared!");
#endif // FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    a_parts.push_back((*add_iter)->mesh_meta_data_ordinal());
  }
  bool quick_verify_check = true;

  for ( AddIterator ia = begin_add_parts; ia != end_add_parts ; ++ia ) {
    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ia, ent_rank, undef_rank);
    const PartVector& supersets = (*ia)->supersets();
    for(PartVector::const_iterator s_iter=supersets.begin(), s_end=supersets.end();
        s_iter!=s_end; ++s_iter) {
      a_parts.push_back((*s_iter)->mesh_meta_data_ordinal());
    }
  }

  order(a_parts);

  OrdinalVector::const_iterator a_parts_begin = a_parts.begin(),
                                a_parts_end   = a_parts.end();
  OrdinalVector r_parts ;

  for ( RemoveIterator ir = begin_remove_parts; ir != end_remove_parts ; ++ir ) {

    // The following guards should be in the public interface to
    // changing parts.  However, internal mechanisms such as changing
    // ownership calls this function to add or remove an entity from
    // the three special parts.  Without refactoring, these guards
    // cannot be put in place.
    /*
    ThrowErrorMsgIf( m_mesh_meta_data.universal_part() == **ir,
                     "Cannot remove entity from universal part" );
    ThrowErrorMsgIf( m_mesh_meta_data.locally_owned_part() == **ir,
                     "Cannot remove entity from locally owned part" );
    ThrowErrorMsgIf( m_mesh_meta_data.globally_shared_part() == **ir,
                     "Cannot remove entity from globally shared part" );
    */

    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ir, ent_rank, undef_rank);

    if ( ! contains_ordinal( a_parts_begin, a_parts_end , (*ir)->mesh_meta_data_ordinal() ) ) {
      r_parts.push_back( (*ir)->mesh_meta_data_ordinal() );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( bucket(entity).member ( **cur_part ) )
          r_parts.push_back ( (*cur_part)->mesh_meta_data_ordinal() );
    }
  }

  order(r_parts);

  // If it looks like we have a problem, run the full check and we should
  // expect to see an exception thrown; otherwise, only do the full check in
  // debug mode because it incurs significant overhead.
  if ( ! quick_verify_check ) {
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
    ThrowRequireMsg(false, "Expected throw from verify methods above.");
  }
  else {
#ifndef NDEBUG
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
#endif
  }

  internal_change_entity_parts( entity , a_parts , r_parts , always_propagate_internal_changes );
}

//  The 'add_parts' and 'remove_parts' are complete and disjoint.
//  Changes need to have parallel resolution during
//  modification_end.

template <typename PartT>
inline
void BulkData::internal_change_entity_parts(
  Entity entity ,
  const std::vector<PartT> & add_parts ,
  const std::vector<PartT> & remove_parts,
  bool always_propagate_internal_changes )
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity);
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "add_parts: " << add_parts);
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "remove_parts: " << remove_parts);

  Bucket * const k_old = bucket_ptr( entity );

  if ( k_old && k_old->member_all( add_parts ) &&
              ! k_old->member_any( remove_parts ) ) {
    // Is already a member of all add_parts,
    // is not a member of any remove_parts,
    // thus nothing to do.
    return ;
  }

  std::vector<PartT> parts_removed ;

  OrdinalVector parts_total ; // The final part list

  //--------------------------------

  if ( k_old ) {
    // Keep any of the existing bucket's parts
    // that are not a remove part.
    // This will include the 'intersection' parts.
    //
    // These parts are properly ordered and unique.

    const std::pair<const unsigned *, const unsigned*>
      bucket_parts = k_old->superset_part_ordinals();

    const unsigned * parts_begin = bucket_parts.first;
    const unsigned * parts_end   = bucket_parts.second;

    const unsigned num_bucket_parts = parts_end - parts_begin;
    parts_total.reserve( num_bucket_parts + add_parts.size() );
    parts_total.insert( parts_total.begin(), parts_begin , parts_end);

    if ( !remove_parts.empty() ) {
      parts_removed.reserve(remove_parts.size());
      impl::filter_out( parts_total , remove_parts , parts_removed );
    }
  }
  else {
    parts_total.reserve(add_parts.size());
  }

  if ( !add_parts.empty() ) {
    impl::merge_in( parts_total , add_parts );
  }

  if ( parts_total.empty() ) {
    // Always a member of the universal part.
    const unsigned univ_ord =
      m_mesh_meta_data.universal_part().mesh_meta_data_ordinal();
    parts_total.push_back( univ_ord );
  }

  EntityRank e_rank = entity_rank(entity);

  //--------------------------------
  // Move the entity to the new partition.
  stk::mesh::impl::Partition *partition =
          m_bucket_repository.get_or_create_partition(e_rank, parts_total);

  if (k_old) {
    k_old->getPartition()->move_to(entity, *partition);
  }
  else {
    partition->add(entity);
  }

  ////
  //// SHOULD WE FIND A WAY TO PUT THE REST OF THIS IN Partition::move_to(..)?
  ////

  // Propagate part changes through the entity's relations.
  //(Only propagate part changes for parts which have a primary-entity-rank that matches
  // the entity's rank. Other parts don't get induced...)

  std::vector<PartT> rank_parts_removed;
  for(typename std::vector<PartT>::const_iterator pr=parts_removed.begin(), prend=parts_removed.end(); pr!=prend; ++pr) {
    if (impl::get_part(*pr, m_mesh_meta_data).primary_entity_rank() == e_rank) {
      rank_parts_removed.push_back(*pr);
    }
  }

  if (always_propagate_internal_changes ||
      !rank_parts_removed.empty() ) {
    internal_propagate_part_changes( entity , rank_parts_removed );
  }
}

//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

template <typename PartT>
inline
void BulkData::internal_propagate_part_changes(
  Entity entity ,
  const std::vector<PartT> & removed )
{
  TraceIfWatching("stk::mesh::BulkData::internal_propagate_part_changes",
      LOG_ENTITY,
      entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity);
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "Removed: " << removed);

  m_check_invalid_rels = false;

  const unsigned erank = entity_rank(entity);
  const EntityRank end_rank = m_mesh_meta_data.entity_rank_count();

  OrdinalVector to_del , to_add , empty ;

  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    size_t num_rels = num_connectivity(entity, irank);
    Entity const *rel_entities = begin(entity, irank);
    ConnectivityOrdinal const *rel_ordinals = begin_ordinals(entity, irank);
    for (size_t j = 0; j < num_rels; ++j)
    {
      const unsigned rel_ident = rel_ordinals[j];

      if ( irank < erank ) { // a 'to' entity

        Entity e_to = rel_entities[j];

        if (e_to == Entity::InvalidEntity)
        {
          continue;
        }

        to_del.clear();
        to_add.clear();
        empty.clear();

        // Induce part membership from this relationship to
        // pick up any additions.
        induced_part_membership(*this, entity, empty, irank, rel_ident, to_add );

        if ( ! removed.empty() ) {
          // Something was removed from the 'from' entity,
          // deduce what may have to be removed from the 'to' entity.

          // Deduce parts for 'e_to' from all upward relations.
          // Any non-parallel part that I removed that is not deduced for
          // 'e_to' must be removed from 'e_to'

          EntityRank e_to_rank = entity_rank(e_to);

          for (EntityRank to_rel_rank_i = e_to_rank + 1; to_rel_rank_i < end_rank; ++to_rel_rank_i)
          {
            size_t num_back_rels = num_connectivity(e_to, to_rel_rank_i);
            Entity const* back_rel_entities = begin(e_to, to_rel_rank_i);
            ConnectivityOrdinal const *rel_ords = begin_ordinals(e_to, to_rel_rank_i);
            for (size_t k = 0; k < num_back_rels; ++k)
            {
              if (entity != back_rel_entities[k])  // Already did this entity
              {
                // Relation from to_rel->entity() to e_to
                induced_part_membership(*this, back_rel_entities[k], empty,
                                        e_to_rank, rel_ords[k], to_add );
              }
            }
          }

          OrdinalVector::const_iterator to_add_begin = to_add.begin(),
            to_add_end   = to_add.end();

          for ( typename std::vector<PartT>::const_iterator
              k = removed.begin() ; k != removed.end() ; ++k ) {
            if ( ! contains_ordinal( to_add_begin, to_add_end , impl::get_ordinal(*k) ) ) {
              induced_part_membership( impl::get_part(*k, m_mesh_meta_data), erank, irank, rel_ident, to_del );
            }
          }
        }

        if ( parallel_size() < 2 || entity_comm_sharing(entity_key(e_to)).empty() ) {
          // Entirely local, ok to remove memberships now
          internal_change_entity_parts( e_to , to_add , to_del );
        }
        else {
          // Shared, do not remove memberships now.
          // Wait until modification_end.
          internal_change_entity_parts( e_to , to_add , empty );
        }
      }
    }
  }
  m_check_invalid_rels = true;
}

// TODO Change the methods below to requirements (private, const invariant checkers)

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 3) Part that does not match the entity rank.

template <typename PartT>
inline
void BulkData::internal_verify_change_parts( const MetaData   & meta ,
                                             const Entity entity ,
                                             const std::vector<PartT> & parts ) const
{
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  const EntityRank undef_rank  = InvalidEntityRank;
  const EntityRank erank = entity_rank(entity);

  bool ok = true ;
  std::ostringstream msg ;

  for ( typename std::vector<PartT>::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part & p = impl::get_part(*i, m_mesh_meta_data);
    const unsigned part_rank = p.primary_entity_rank();

    bool intersection_ok, rel_target_ok, rank_ok;
    internal_basic_part_check(&p, erank, undef_rank, intersection_ok, rel_target_ok, rank_ok);

    if ( !intersection_ok || !rel_target_ok || !rank_ok ) {
      if ( ok ) {
        ok = false ;
        msg << "change parts for entity " << identifier(entity);
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }

      msg << p.name() << "[" ;
      if ( part_rank < rank_names.size() ) {
        msg << rank_names[ part_rank ];
      }
      else {
        msg << part_rank ;
      }
      msg << "] " ;
      if ( !intersection_ok ) { msg << "is_intersection " ; }
      if ( !rel_target_ok )   { msg << "is_relation_target " ; }
      if ( !rank_ok )         { msg << "is_bad_rank " ; }
    }
  }

  ThrowErrorMsgIf( !ok, msg.str() << "}" );
}

} //namespace mesh
} //namespace stk;

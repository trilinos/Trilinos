#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

std::ostream &
print_relation( std::ostream & s , relation_attr_type attr )
{
  s << relation_entity_type( attr ); //this used to be entity_type_name(...)
  s << "[" ;
  s << relation_kind( attr );
  s << "." ;
  s << relation_identifier( attr );
  s << "]" ;
  return s ;
}

std::ostream &
print_relation( std::ostream & s , const MetaData & meta ,
                relation_attr_type attr , EntityKey key )
{
  print_relation( s , attr );
  print_entity_key( s , meta , key );

  return s ;
}

std::ostream &
operator << ( std::ostream & s , const Relation & con )
{
  Entity * const e = con.entity();

  print_relation( s , con.attribute() );

  if ( e ) {
    print_entity_key( s , e->bucket().mesh().mesh_meta_data(), e->key() );
  }
  else {
    s << "?" ;
  }

  return s ;
}

//----------------------------------------------------------------------

Relation::Relation( Entity & entity ,
                    unsigned identifier ,
                    unsigned kind )
  : m_attr( relation_attr( entity.entity_type() , identifier , kind ) ),
    m_entity( & entity )
{}

Relation::Relation( relation_attr_type attr , Entity & entity )
  : m_attr( attr ), m_entity( & entity )
{
  if ( relation_entity_type( attr ) != entity.entity_type() ) {
    std::ostringstream msg ;
    msg << "stk::mesh::Relation::Relation( "  ;
    print_relation( msg , attr );
    msg << " , " ;
    print_entity_key( msg , entity.bucket().mesh().mesh_meta_data() , entity.key() );
    msg << " ) INCOMPATIBLE ARGUMENTS" ;
    throw std::invalid_argument( msg.str() );
  }
}

bool Relation::operator < ( const Relation & r ) const
{
  bool result = false;

  if ( m_attr != r.m_attr ) {
    result = m_attr < r.m_attr ;
  }
  else {
    const EntityKey lhs = m_entity   ? m_entity->key()   : EntityKey() ;
    const EntityKey rhs = r.m_entity ? r.m_entity->key() : EntityKey() ;
    result = lhs < rhs ;
  }
  return result ;
}

//----------------------------------------------------------------------

namespace {

struct LessRelation {

  inline
  bool operator()( const Relation & lhs , const Relation & rhs ) const
    { return lhs.operator < ( rhs ); }

  inline
  bool operator()( const Relation & lhs , const relation_attr_type rhs ) const
    { return lhs.attribute() < rhs ; }
};
}

PairIterRelation Entity::relations( unsigned type , unsigned kind ) const
{
  const unsigned et_next = type + 1 ;
  const relation_attr_type lo_attr = relation_attr( type ,    0 , kind ),
                           hi_attr = relation_attr( et_next , 0 , kind );

  std::vector<Relation>::const_iterator i = m_relation.begin();
  std::vector<Relation>::const_iterator e = m_relation.end();

  i = std::lower_bound( i , e , lo_attr , LessRelation() );
  e = std::lower_bound( i , e , hi_attr , LessRelation() );

  return PairIterRelation( i , e );
}


//----------------------------------------------------------------------

/** \brief  Query if a member entity of the given entity type 
 *          has an induced membership.
 */
bool membership_is_induced( const Part & part , unsigned entity_type )
{
  const MetaData & meta = part.mesh_meta_data();

  const bool induced_by_type =
     entity_type < part.primary_entity_type() &&
                   part.primary_entity_type() < meta.entity_type_count() ;
 
  const bool induced_by_stencil =
    ! part.relations().empty() &&
      part.relations().begin()->m_target == & part ;
 
  return induced_by_type || induced_by_stencil ;
}
 
//----------------------------------------------------------------------

void induced_part_membership( Part & part , 
                              unsigned entity_type_from , 
                              unsigned entity_type_to ,
                              unsigned relation_identifier ,
                              unsigned relation_kind ,
                              PartVector & induced_parts )
{
  if ( entity_type_to < entity_type_from &&
       part.primary_entity_type() == entity_type_from ) {
 
    // Direct relationship:
 
    insert( induced_parts , part );
 
    // Stencil relationship where 'part' is the root:
    // The 'target' should not have subsets or supersets.
 
    const std::vector<PartRelation> & part_rel = part.relations();
 
    for ( std::vector<PartRelation>::const_iterator
          j = part_rel.begin() ; j != part_rel.end() ; ++j ) {
 
      if ( & part == j->m_root &&
           0 <= (* j->m_function)( entity_type_from , entity_type_to ,
                                   relation_identifier , relation_kind ) ) {
        insert( induced_parts , * j->m_target );
      }
    }
  }
}
 
//----------------------------------------------------------------------
//  What are this entity's part memberships that can be deduced from
//  this entity's relationship.  Can only trust 'entity_from' to be
//  accurate if it is owned by the local process.

void induced_part_membership( const Entity     & entity_from ,
                              const PartVector & omit ,
                                    unsigned     entity_type_to ,
                                    unsigned     relation_identifier ,
                                    unsigned     relation_kind ,
                                    PartVector & entity_to_parts )
{
  const Bucket   & bucket_from    = entity_from.bucket();
  const BulkData & mesh           = bucket_from.mesh();
  const unsigned local_proc_rank  = mesh.parallel_rank();
  const unsigned entity_type_from = entity_from.entity_type();

  if ( entity_type_to < entity_type_from &&
       local_proc_rank == entity_from.owner_rank() ) {
    const MetaData   & meta        = mesh.mesh_meta_data();
    const PartVector & all_parts   = meta.get_parts();

    const std::pair<const unsigned *, const unsigned *>
      bucket_superset_ordinals = bucket_from.superset_part_ordinals();

    // Contributions of the 'from' entity:

    for ( const unsigned * i = bucket_superset_ordinals.first ;
                           i < bucket_superset_ordinals.second ; ++i ) {

      Part & part = * all_parts[*i] ;

      if ( ! contain( omit , part ) ) {
        induced_part_membership( part,
                                 entity_type_from ,
                                 entity_type_to ,
                                 relation_identifier ,
                                 relation_kind , entity_to_parts );
      }
    }
  }
}

//----------------------------------------------------------------------

void induced_part_membership( const Entity     & entity ,
                              const PartVector & omit ,
                                    PartVector & induced )
{
  for ( PairIterRelation 
        rel = entity.relations() ; ! rel.empty() ; ++rel ) {  

    induced_part_membership( * rel->entity() , omit ,
                             entity.entity_type() ,
                             rel->identifier() ,
                             rel->kind() ,
                             induced );
  }
}

//----------------------------------------------------------------------

void set_field_relations( Entity & e_from ,
                          Entity & e_to ,
                          const unsigned ident ,
                          const unsigned kind )
{
  const std::vector<FieldRelation> & field_rels =
    e_from.bucket().mesh().mesh_meta_data().get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator 
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      void * const src = field_data( * fr.m_target , e_to );

      const size_t number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const size_t offset =
         (*fr.m_function)( e_from.entity_type() , 
                           e_to.entity_type() , ident , kind );

      if ( offset < number ) {
        ptr[ offset ] = src ;
      }
    }
  }
}

//----------------------------------------------------------------------

namespace {

void assert_valid_relation( const char method[] ,
                            const BulkData & mesh ,
                            const Entity   & e_from ,
                            const Entity   & e_to )
{
  const bool error_mesh_from = & mesh != & e_from.bucket().mesh();
  const bool error_mesh_to   = & mesh != & e_to.bucket().mesh();
  const bool error_type      = e_from.entity_type() <= e_to.entity_type();

  if ( error_mesh_from || error_mesh_to || error_type ) {
    std::ostringstream msg ;
    msg << method << "( from " ;
    print_entity_key( msg , mesh.mesh_meta_data(), e_from.key() );
    if ( error_mesh_from ) {
      msg << " NOT MEMBER OF THIS MESH" ;
    }
    msg << " , to " ;
    print_entity_key( msg , mesh.mesh_meta_data(), e_to.key() );
    if ( error_mesh_to ) {
      msg << " NOT MEMBER OF THIS MESH" ;
    }
    msg << " ) FAILED" ;
    if ( error_type ) {
      msg << " A relation must be from higher to lower ranking entity" ;
    }
    throw std::runtime_error( msg.str() );
  }
}

}

void BulkData::declare_relation( Entity & e_from ,
                                 Entity & e_to ,
                                 const unsigned local_id ,
                                 const unsigned kind )
{
  static const char method[] = "stk::mesh::BulkData::declare_relation" ;

  assert_ok_to_modify( method );

  assert_valid_relation( method , *this , e_from , e_to );

  bool error_internal = false ;

  const Relation forward( e_to , local_id , kind );

  const std::vector<Relation>::iterator fe = e_from.m_relation.end();
        std::vector<Relation>::iterator fi = e_from.m_relation.begin();

  fi = std::lower_bound( fi , fe , forward , LessRelation() );

  if ( fe == fi || forward.attribute() != fi->attribute() ) {

    // A new relation and its converse

    const Relation converse( e_from , local_id , kind );
    const std::vector<Relation>::iterator ce = e_to.m_relation.end();
          std::vector<Relation>::iterator ci = e_to.m_relation.begin();

    ci = std::lower_bound( ci , ce , converse , LessRelation() );

    if ( ce == ci || converse != *ci ) {
      fi = e_from.m_relation.insert( fi , forward );
      ci = e_to  .m_relation.insert( ci , converse );

      e_from.m_sync_count = m_sync_count ;
      e_to  .m_sync_count = m_sync_count ;

      PartVector add , empty ;

      // Deduce and set new part memberships:

      induced_part_membership( e_from, empty,
                               e_to.entity_type(), local_id, kind, add );

      internal_change_entity_parts( e_to , add , empty );

      set_field_relations( e_from , e_to , local_id , kind );
    }
    else {
      error_internal = true ;
    }
  }

  if ( error_internal || fi->entity() != & e_to ) {
    std::ostringstream msg ;
    msg << method << "( from " ;
    print_entity_key( msg , m_mesh_meta_data, e_from.key() );
    msg << " , to " ;
    print_entity_key( msg , m_mesh_meta_data, e_to.key() );
    msg << " , id " << local_id ;
    msg << " , kind " << kind ;
    msg << " ) FAILED" ;
    if ( fi->entity() != & e_to ) {
      msg << " Relation already exists to " ;
      print_entity_key( msg , m_mesh_meta_data, fi->entity()->key() );
    }
    if ( error_internal ) {
      msg << " Internal error - converse relation already exists" ;
    }
    throw std::runtime_error( msg.str() );
  }
}

void BulkData::declare_relation( Entity & entity ,
                                 const std::vector<Relation> & rel )
{
  static const char method[] = "stk::mesh::BulkData::declare_relation" ;

  assert_ok_to_modify( method );

  const unsigned etype = entity.entity_type();

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity & e = * i->entity();
    const unsigned k = i->kind();
    const unsigned n = i->identifier();
    if ( e.entity_type() < etype ) {
      declare_relation( entity , e , n , k );
    }
    else if ( etype < e.entity_type() ) {
      declare_relation( e , entity , n , k );
    }
    else {
      throw std::runtime_error( std::string("stk::mesh::BulkData::declare_relation FAILED: given entities of the same entity type") ); 
    }
  }
}

//----------------------------------------------------------------------

namespace {

void clear_field_relations( Entity & e_from ,
                            const unsigned type ,
                            const unsigned ident ,
                            const unsigned kind )
                          
{
  const std::vector<FieldRelation> & field_rels =
    e_from.bucket().mesh().mesh_meta_data().get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator 
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      const size_t number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const size_t offset =
        (*fr.m_function)( e_from.entity_type() , type , ident , kind );

      if ( offset < number ) {
        ptr[ offset ] = NULL ;
      }
    }
  }
}

}

void BulkData::destroy_relation( Entity & e_from ,
                                 Entity & e_to , unsigned kind )
{
  static const char method[]= "stk::mesh::BulkData::destroy_relation" ;

  assert_ok_to_modify( method );

  assert_valid_relation( method , *this , e_from , e_to );

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  PartVector del , keep ;

  for ( std::vector<Relation>::iterator
        i = e_to.m_relation.end() ; i != e_to.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e_from && i->kind() == kind ) {
      i = e_to.m_relation.erase( i );
    }
    else if ( e_to.entity_type() < i->entity_type() ) {
      induced_part_membership( * i->entity(), del, e_to.entity_type(),
                               i->identifier(), i->kind(), keep );
    }
  }

  for ( std::vector<Relation>::iterator
        i = e_from.m_relation.end() ; i != e_from.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e_to && i->kind() == kind ) {

      induced_part_membership( e_from, keep, e_to.entity_type(),
                               i->identifier(), i->kind(), del );

      clear_field_relations( e_from , e_to.entity_type() ,
                             i->identifier() , i->kind() );

      i = e_from.m_relation.erase( i );
    }
  }

  //------------------------------
  // 'keep' contains the parts deduced from kept relations
  // 'del'  contains the parts deduced from deleted relations
  //        that are not in 'keep'
  // Only remove these part memberships the entity is not shared.
  // If the entity is shared then wait until modificaton_end_synchronize.
  //------------------------------

  if ( ! del.empty() && e_to.sharing().empty() ) {
    PartVector add ;
    internal_change_entity_parts( e_to , add , del );
  }
}

//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

void BulkData::internal_propagate_part_changes(
  Entity           & entity ,
  const PartVector & removed )
{
  const unsigned etype = entity.entity_type();

  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const unsigned rel_type  = rel->entity_type();
    const unsigned rel_ident = rel->identifier();
    const unsigned rel_kind  = rel->kind();

    if ( rel_type < etype ) { // a 'to' entity

      Entity & e_to = * rel->entity();

      PartVector to_del , to_add , empty ;

      // Induce part membership from this relationship to
      // pick up any additions.
      induced_part_membership( entity, empty,
                               rel_type, rel_ident, rel_kind, to_add );

      if ( ! removed.empty() ) {
        // Something was removed from the 'from' entity,
        // deduce what may have to be removed from the 'to' entity.

        // Deduce parts for 'e_to' from all upward relations.
        // Any non-parallel part that I removed that is not deduced for
        // 'e_to' must be removed from 'e_to'

        for ( PairIterRelation
              to_rel = e_to.relations(); ! to_rel.empty() ; ++to_rel ) {
          if ( e_to.entity_type() < to_rel->entity_type() &&
               & entity != to_rel->entity() /* Already did this entity */ ) {
            // Relation from to_rel->entity() to e_to
            induced_part_membership( * to_rel->entity(), empty,
                                     e_to.entity_type(),
                                     to_rel->identifier(),
                                     to_rel->kind(), to_add );
          }
        }

        for ( PartVector::const_iterator
              j = removed.begin() ; j != removed.end() ; ++j ) {
          if ( ! contain( to_add , **j ) ) {
            induced_part_membership( **j, etype, rel_type, rel_ident, rel_kind,
                                     to_del );
          }
        }
      }

      if ( e_to.sharing().empty() ) {
        // Entirely local, ok to remove memberships now
        internal_change_entity_parts( e_to , to_add , to_del );
      }
      else {
        // Shared, do not remove memberships now.
        // Wait until modification_end.
        internal_change_entity_parts( e_to , to_add , empty );
      }

      set_field_relations( entity, e_to, rel_ident , rel_kind );
    }
    else if ( etype < rel_type ) { // a 'from' entity
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel_ident, rel_kind );
    }
  }
}

void BulkData::internal_propagate_relocation( Entity & entity )
{
  const unsigned etype = entity.entity_type();
  PairIterRelation rel = entity.relations();

  for ( ; ! rel.empty() ; ++rel ) {
    const unsigned rel_type = rel->entity_type();
    if ( rel_type < etype ) {
      Entity & e_to = * rel->entity();

      set_field_relations( entity, e_to, rel->identifier(), rel->kind() );
    }
    else if ( etype < rel_type ) {
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel->identifier(), rel->kind() );
    }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


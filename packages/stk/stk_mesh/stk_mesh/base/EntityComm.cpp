/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iterator>
#include <stdexcept>
#include <sstream>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/EntityComm.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------------

bool in_shared( const Entity & entity )
{
  PairIterEntityComm ec = entity.comm();
  return ! ec.empty() && ec.front().ghost_id == 0 ;
}

bool in_shared( const Entity & entity , unsigned proc )
{
  for ( PairIterEntityComm ec = entity.comm();
        ! ec.empty() && ec->ghost_id == 0 ; ++ec ) {
    if ( proc == ec->proc ) {
      return true ;
    }
  }
  return false ;
}

bool in_receive_ghost( const Entity & entity )
{
  // Ghost communication with owner.
  PairIterEntityComm ec = entity.comm();
  return ! ec.empty() && ec.front().ghost_id != 0 &&
                         ec.front().proc == entity.owner_rank();
}

bool in_receive_ghost( const Ghosting & ghost , const Entity & entity )
{
  return in_ghost( ghost , entity , entity.owner_rank() );
}

bool in_send_ghost( const Entity & entity )
{
  // Ghost communication with non-owner.
  PairIterEntityComm ec = entity.comm();
  return ! ec.empty() && ec.back().ghost_id != 0 &&
                         ec.back().proc != entity.owner_rank();
}

bool in_send_ghost( const Entity & entity , unsigned proc )
{
  for ( PairIterEntityComm ec = entity.comm(); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id != 0 &&
         ec->proc   != entity.owner_rank() &&
         ec->proc   == proc ) {
      return true ;
    }
  }
  return false ;
}

bool in_ghost( const Ghosting & ghost , const Entity & entity , unsigned p )
{
  // Ghost communication from owner.
  EntityCommInfo tmp( ghost.ordinal() , p );

  std::vector<EntityCommInfo>::const_iterator i =
    std::lower_bound( entity.comm().begin() , entity.comm().end() , tmp );

  return i != entity.comm().end() && tmp == *i ;
}

/** \brief  Is in owned closure of the given process,
 *          typically the local process.
 */
bool in_owned_closure( const Entity & entity , unsigned proc )
{
  bool result = entity.owner_rank() == proc ;

  if ( ! result ) {
    const unsigned erank = entity.entity_rank();

    for ( PairIterRelation
          rel = entity.relations(); ! result && ! rel.empty() ; ++rel ) {
      result =  erank < rel->entity_rank() &&
                in_owned_closure( * rel->entity(), proc);
    }
  }

  return result ;
}

void comm_procs( const Entity & entity , std::vector<unsigned> & procs )
{
  procs.clear();
  for ( PairIterEntityComm ec = entity.comm(); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<unsigned>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}

void comm_procs( const Ghosting & ghost ,
                 const Entity & entity , std::vector<unsigned> & procs )
{
  procs.clear();
  for ( PairIterEntityComm ec = entity.comm(); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id == ghost.ordinal() ) {
      procs.push_back( ec->proc );
    }
  }
}


//----------------------------------------------------------------------------

void pack_entity_info( CommBuffer & buf , const Entity & entity )
{
  const EntityKey & key   = entity.key();
  const unsigned    owner = entity.owner_rank();
  const std::pair<const unsigned *, const unsigned *>
    part_ordinals = entity.bucket().superset_part_ordinals();
  const PairIterRelation relations = entity.relations();

  const unsigned nparts = part_ordinals.second - part_ordinals.first ;
  const unsigned nrel   = relations.size();

  buf.pack<EntityKey>( key );
  buf.pack<unsigned>( owner );
  buf.pack<unsigned>( nparts );
  buf.pack<unsigned>( part_ordinals.first , nparts );
  buf.pack<unsigned>( nrel );

  for ( unsigned i = 0 ; i < nrel ; ++i ) {
    buf.pack<EntityKey>( relations[i].entity()->key() );
    buf.pack<Relation::raw_attr_type>( relations[i].attribute() );
  }
}

void unpack_entity_info(
  CommBuffer       & buf,
  const BulkData   & mesh ,
  EntityKey        & key ,
  unsigned         & owner ,
  PartVector       & parts ,
  std::vector<Relation> & relations )
{
  unsigned nparts = 0 ;
  unsigned nrel = 0 ;

  buf.unpack<EntityKey>( key );
  buf.unpack<unsigned>( owner );
  buf.unpack<unsigned>( nparts );

  parts.resize( nparts );

  for ( unsigned i = 0 ; i < nparts ; ++i ) {
    unsigned part_ordinal = ~0u ;
    buf.unpack<unsigned>( part_ordinal );
    parts[i] = & mesh.mesh_meta_data().get_part( part_ordinal );
  }

  buf.unpack( nrel );

  relations.clear();
  relations.reserve( nrel );

  for ( unsigned i = 0 ; i < nrel ; ++i ) {
    EntityKey rel_key ;
    Relation::raw_attr_type rel_attr = 0 ;
    buf.unpack<EntityKey>( rel_key );
    buf.unpack<Relation::raw_attr_type>( rel_attr );
    Entity * const entity =
      mesh.get_entity( entity_rank(rel_key), entity_id(rel_key) );
    if ( entity && EntityLogDeleted != entity->log_query() ) {
      Relation rel( rel_attr , * entity );
      relations.push_back( rel );
    }
  }
}


//----------------------------------------------------------------------

void pack_field_values( CommBuffer & buf , Entity & entity )
{
  const Bucket   & bucket = entity.bucket();
  const BulkData & mesh   = bucket.mesh();
  const MetaData & mesh_meta_data = mesh.mesh_meta_data();

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  for ( std::vector< FieldBase * >::const_iterator
        i = fields.begin() ; i != fields.end() ; ++i ) {

    const FieldBase & f = **i ;

    if ( f.data_traits().is_pod ) {
      const unsigned size = field_data_size( f , bucket );

      buf.pack<unsigned>( size );

      if ( size ) {
        unsigned char * const ptr =
          reinterpret_cast<unsigned char *>( field_data( f , entity ) );
        buf.pack<unsigned char>( ptr , size );
      }
    }
  }
}

bool unpack_field_values(
  CommBuffer & buf , Entity & entity , std::ostream & error_msg )
{
  const Bucket   & bucket = entity.bucket();
  const BulkData & mesh   = bucket.mesh();
  const MetaData & mesh_meta_data = mesh.mesh_meta_data();

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  bool ok = true ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;

    if ( f.data_traits().is_pod ) {

      const unsigned size = field_data_size( f , bucket );
      unsigned recv_data_size = 0 ;
      buf.unpack<unsigned>( recv_data_size );

      if ( size != recv_data_size ) {
        if ( ok ) {
          ok = false ;
          print_entity_key( error_msg , mesh_meta_data , entity.key() );
        }
        error_msg << " " << f.name();
        error_msg << " " << size ;
        error_msg << " != " << recv_data_size ;
        buf.skip<unsigned char>( recv_data_size );
      }
      else if ( size ) { // Non-zero and equal
        unsigned char * ptr =
          reinterpret_cast<unsigned char *>( field_data( f , entity ) );
        buf.unpack<unsigned char>( ptr , size );
      }
    }
  }

  return ok ;
}

//----------------------------------------------------------------------

} // namespace mesh
}


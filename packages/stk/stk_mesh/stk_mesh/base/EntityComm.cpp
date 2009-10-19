#include <iterator>
#include <stdexcept>
#include <sstream>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/EntityComm.hpp>

namespace stk {
namespace mesh {

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
    buf.pack<relation_attr_type>( relations[i].attribute() );
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
    relation_attr_type rel_attr = 0 ;
    buf.unpack<EntityKey>( rel_key );
    buf.unpack<relation_attr_type>( rel_attr );
    Entity * const entity = mesh.get_entity( entity_type(rel_key), entity_id(rel_key) );
    if ( entity ) {
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

    const unsigned size = field_data_size( f , bucket );

    buf.pack<unsigned>( size );

    if ( size ) {
      unsigned char * const ptr =
        reinterpret_cast<unsigned char *>( field_data( f , entity ) );
      buf.pack<unsigned char>( ptr , size );
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

  return ok ;
}

//----------------------------------------------------------------------

} // namespace mesh
}


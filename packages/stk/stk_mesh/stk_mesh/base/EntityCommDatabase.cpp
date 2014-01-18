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
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------------

void pack_entity_info(const BulkData& mesh, CommBuffer & buf , const Entity entity )
{
  const EntityKey & key   = mesh.entity_key(entity);
  const unsigned    owner = mesh.parallel_owner_rank(entity);
  const std::pair<const unsigned *, const unsigned *>
    part_ordinals = mesh.bucket(entity).superset_part_ordinals();

  const unsigned nparts = part_ordinals.second - part_ordinals.first ;
  const unsigned tot_rel = mesh.count_relations(entity);
  Bucket& bucket = mesh.bucket(entity);
  unsigned ebo   = mesh.bucket_ordinal(entity);

  buf.pack<EntityKey>( key );
  buf.pack<unsigned>( owner );
  buf.pack<unsigned>( nparts );
  buf.pack<unsigned>( part_ordinals.first , nparts );
  buf.pack<unsigned>( tot_rel );

  ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity " << entity.local_offset());
  const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();

  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    const unsigned nrel = bucket.num_connectivity(ebo, irank);
    if (nrel > 0) {
      Entity const *rel_entities = bucket.begin(ebo, irank);
      ConnectivityOrdinal const *rel_ordinals = bucket.begin_ordinals(ebo, irank);
      Permutation const *rel_permutations = bucket.begin_permutations(ebo, irank);

      for ( unsigned i = 0 ; i < nrel ; ++i ) {
        if (mesh.is_valid(rel_entities[i])) {
          buf.pack<EntityKey>( mesh.entity_key(rel_entities[i]) );
          buf.pack<unsigned>( rel_ordinals[i] );
          if (bucket.has_permutation(irank)) {
            buf.pack<unsigned>( rel_permutations[i] );
          } else {
            buf.pack<unsigned>(0u);
          }
        } else { // relation to invalid entity (FIXED CONNECTIVITY CASE)
          // TODO:  Consider not communicating relations to invalid entities...
          buf.pack<EntityKey>( EntityKey() ); // invalid EntityKey
          buf.pack<unsigned>( rel_ordinals[i] );
          buf.pack<unsigned>(0u); // permutation
        }
      }
    }
  }
}

void unpack_entity_info(
  CommBuffer       & buf,
  const BulkData   & mesh ,
  EntityKey        & key ,
  int              & owner ,
  PartVector       & parts ,
  std::vector<Relation> & relations )
{
  unsigned nparts = 0 ;
  unsigned nrel = 0 ;

  buf.unpack<EntityKey>( key );
  buf.unpack<int>( owner );
  buf.unpack<unsigned>( nparts );

  parts.resize( nparts );

  for ( unsigned i = 0 ; i < nparts ; ++i ) {
    unsigned part_ordinal = ~0u ;
    buf.unpack<unsigned>( part_ordinal );
    parts[i] = & MetaData::get(mesh).get_part( part_ordinal );
  }

  buf.unpack( nrel );

  relations.clear();
  relations.reserve( nrel );

  for ( unsigned i = 0 ; i < nrel ; ++i ) {
    EntityKey rel_key ;
    unsigned rel_id = 0 ;
    unsigned rel_attr = 0 ;
    buf.unpack<EntityKey>( rel_key );
    buf.unpack<unsigned>( rel_id );
    buf.unpack<unsigned>( rel_attr );
    if (rel_key == EntityKey()) {
      continue;
    }
    Entity const entity =
      mesh.get_entity( rel_key.rank(), rel_key.id() );
    if ( mesh.is_valid(entity) ) {
      Relation rel(mesh, entity, rel_id );
      rel.set_attribute(rel_attr);
      relations.push_back( rel );
    }
  }
}


//----------------------------------------------------------------------

void pack_field_values(const BulkData& mesh, CommBuffer & buf , Entity entity )
{
  if (!mesh.is_field_updating_active()) {
    return;
  }

  const Bucket   & bucket = mesh.bucket(entity);
  const MetaData & mesh_meta_data = MetaData::get(mesh);

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  for ( std::vector< FieldBase * >::const_iterator
        i = fields.begin() ; i != fields.end() ; ++i ) {

    const FieldBase & f = **i ;

    if(mesh.is_matching_rank(f, bucket)) {


      if ( f.data_traits().is_pod ) {
	const unsigned size = mesh.field_data_size_per_entity( f, bucket );

	buf.pack<unsigned>( size );

	if ( size ) {
	  unsigned char * const ptr =
	    reinterpret_cast<unsigned char *>( mesh.field_data( f , entity ) );
	  buf.pack<unsigned char>( ptr , size );
	}
      }
    }
  }
}

bool unpack_field_values(const BulkData& mesh,
  CommBuffer & buf , Entity entity , std::ostream & error_msg )
{
  if (!mesh.is_field_updating_active()) {
    return true;
  }

  const Bucket   & bucket = mesh.bucket(entity);
  const MetaData & mesh_meta_data = MetaData::get(mesh);

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  bool ok = true ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;

    if(mesh.is_matching_rank(f, bucket)) {

      if ( f.data_traits().is_pod ) {

	const unsigned size = mesh.field_data_size_per_entity( f, bucket );
	unsigned recv_data_size = 0 ;
	buf.unpack<unsigned>( recv_data_size );

	if ( size != recv_data_size ) {
	  if ( ok ) {
	    ok = false ;
	    error_msg << mesh.identifier(entity);
	  }
	  error_msg << " " << f.name();
	  error_msg << " " << size ;
	  error_msg << " != " << recv_data_size ;
	  buf.skip<unsigned char>( recv_data_size );
	}
	else if ( size ) { // Non-zero and equal
	  unsigned char * ptr =
	    reinterpret_cast<unsigned char *>( mesh.field_data( f , entity ) );
	  buf.unpack<unsigned char>( ptr , size );
	}
      }
    }
  }

  return ok ;
}

//----------------------------------------------------------------------

} // namespace mesh
}


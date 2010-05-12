/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------
#include <stdexcept>
#include <stk_mesh/baseImpl/BucketImpl.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

  
//----------------------------------------------------------------------
namespace {

void * local_malloc( size_t n )
{
  void * const ptr = malloc( n );

  if ( NULL == ptr ) {
    std::ostringstream msg ;
    msg << "stk::mesh::impl::BucketImpl::declare_bucket FAILED malloc( " << n << " )" ;
    throw std::runtime_error( msg.str() );
  }

  return ptr ;
}


void memory_copy( unsigned char * dst , unsigned char * src , unsigned n )
{ memcpy( dst , src , n ); }

void memory_zero( unsigned char * dst , unsigned n )
{ memset( dst , 0 , n ); }

} // namespace

//----------------------------------------------------------------------
BucketImpl::BucketImpl( BulkData        & arg_mesh ,
                unsigned          arg_entity_rank ,
                const unsigned  * arg_key ,
                size_t            arg_alloc_size ,
                size_t            arg_capacity ,
                impl::BucketImpl::DataMap * arg_field_map ,
                Entity         ** arg_entity_array )
: m_mesh( arg_mesh ) ,
  m_entity_rank( arg_entity_rank ) ,
  m_key( arg_key ) ,
  m_alloc_size( arg_alloc_size ) ,
  m_capacity( arg_capacity ) ,
  m_size( 0 ) ,
  m_bucket() ,
  m_field_map( arg_field_map ) ,
  m_entities( arg_entity_array )
{}


//----------------------------------------------------------------------

namespace {

inline unsigned align( size_t nb )
{
  enum { BYTE_ALIGN = 16 };
  const unsigned gap = nb % BYTE_ALIGN ;
  if ( gap ) { nb += BYTE_ALIGN - gap ; }
  return nb ;
}

struct FieldRestrictionLess {
  bool operator()( const FieldBase::Restriction & lhs ,
                   const EntityKey & rhs ) const
    { return lhs.key < rhs ; }
};

const FieldBase::Restriction & empty_field_restriction()
{
  static const FieldBase::Restriction empty ;
  return empty ;
}

const FieldBase::Restriction & dimension( const FieldBase & field ,
                                          unsigned etype ,
                                          const unsigned num_part_ord ,
                                          const unsigned part_ord[] ,
                                          const char * const method )
{
  const FieldBase::Restriction & empty = empty_field_restriction();
  const FieldBase::Restriction * dim = & empty ;

  const std::vector<FieldBase::Restriction> & dim_map = field.restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator iend = dim_map.end();
        std::vector<FieldBase::Restriction>::const_iterator ibeg = dim_map.begin();

  for ( unsigned i = 0 ; i < num_part_ord && iend != ibeg ; ++i ) {

    const EntityKey key = EntityKey(etype,part_ord[i]);

    ibeg = std::lower_bound( ibeg , iend , key , FieldRestrictionLess() );

    if ( iend != ibeg && ibeg->key == key ) {
      if ( dim == & empty ) { dim = & *ibeg ; }

      if ( Compare< MaximumFieldDimension >::
             not_equal( ibeg->stride , dim->stride ) ) {

        Part & p_old = field.mesh_meta_data().get_part( ibeg->ordinal() );
        Part & p_new = field.mesh_meta_data().get_part( dim->ordinal() );

        std::ostringstream msg ;
        msg << method ;
        msg << " FAILED WITH INCOMPATIBLE DIMENSIONS FOR " ;
        msg << field ;
        msg << " Part[" << p_old.name() ;
        msg << "] and Part[" << p_new.name() ;
        msg << "]" ;

        throw std::runtime_error( msg.str() );
      }
    }
  }

  return *dim ;
}

} // namespace

//----------------------------------------------------------------------
// The current 'last' bucket in a family is to be deleted.
// The previous 'last' bucket becomes the new 'last' bucket in the family.

void BucketImpl::destroy_bucket( std::vector<Bucket*> & bucket_set , Bucket * last )
{
  static const char method[] = "stk::mesh::impl::BucketImpl::destroy_bucket" ;

  Bucket * const first = bucket_counter( last->key() ) ? last->m_bucketImpl.m_bucket : last ;

  if ( 0 != last->size() || last != first->m_bucketImpl.m_bucket ) {
    throw std::logic_error(std::string(method));
  }

  std::vector<Bucket*>::iterator ik = lower_bound(bucket_set, last->key());

  if ( ik == bucket_set.end() ) {
    throw std::logic_error(std::string(method));
  }

  if ( last != *ik ) {
    throw std::logic_error(std::string(method));
  }

  ik = bucket_set.erase( ik );

  if ( first != last ) {

    if ( ik == bucket_set.begin() ) {
      throw std::logic_error(std::string(method));
    }

    first->m_bucketImpl.m_bucket = *--ik ;

    if ( 0 == first->m_bucketImpl.m_bucket->size() ) {
      throw std::logic_error(std::string(method));
    }
  }

  BucketImpl::destroy_bucket( last );
}

//----------------------------------------------------------------------
void BucketImpl::destroy_bucket( Bucket * k )
{
  if ( 0 == bucket_counter( k->key() ) ) {
    free( k->m_bucketImpl.m_field_map );
  }

  k->~Bucket();

  free( k );
}


//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.
Bucket *
BucketImpl::declare_nil_bucket( BulkData & mesh , unsigned field_count )
{
  //----------------------------------
  // Field map gives NULL for all field data.

  impl::BucketImpl::DataMap * field_map =
    reinterpret_cast<impl::BucketImpl::DataMap*>(
      local_malloc( sizeof(impl::BucketImpl::DataMap) * ( field_count + 1 )));

  const FieldBase::Restriction & dim = empty_field_restriction();

  for ( unsigned i = 0 ; i < field_count ; ++i ) {
    field_map[ i ].m_base = 0 ;
    field_map[ i ].m_size = 0 ;
    field_map[ i ].m_stride = dim.stride ;
  }
  field_map[ field_count ].m_base   = 0 ;
  field_map[ field_count ].m_size   = 0 ;
  field_map[ field_count ].m_stride = NULL ;

  //----------------------------------
  // Allocation size:  sizeof(Bucket) + key_size * sizeof(unsigned);

  const unsigned alloc_size = align( sizeof(Bucket) ) +
                              align( sizeof(unsigned) * 2 );

  // All fields checked and sized, Ready to allocate

  void * const alloc_ptr = local_malloc( alloc_size );

  unsigned char * ptr = reinterpret_cast<unsigned char *>( alloc_ptr );

  ptr += align( sizeof( Bucket ) );

  unsigned * const new_key = reinterpret_cast<unsigned *>( ptr );

  // Key layout:
  // { part_count + 1 , { part_ordinals } , family_count }

  new_key[0] = 1 ; // part_count + 1
  new_key[1] = 0 ; // family_count

  const unsigned bad_entity_rank = ~0u ;

  Bucket * bucket =
    new( alloc_ptr ) Bucket( mesh , bad_entity_rank , new_key ,
                             alloc_size , 0 , field_map , NULL );

  bucket->m_bucketImpl.m_bucket = bucket ;

  //----------------------------------

  return bucket ;
}


//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.
Bucket *
BucketImpl::declare_bucket( BulkData & mesh ,
                        const unsigned arg_entity_rank ,
                        const unsigned part_count ,
                        const unsigned part_ord[] ,
                        const unsigned bucket_capacity ,
                        const std::vector< FieldBase * > & field_set ,
                              std::vector<Bucket*>       & bucket_set )
{
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  static const char method[] = "stk::mesh::impl::BucketImpl::declare_bucket" ;

  const unsigned max = ~(0u);
  const size_t   num_fields = field_set.size();

  //----------------------------------
  // For performance try not to allocate a temporary.

  unsigned key_tmp_buffer[ KEY_TMP_BUFFER_SIZE ];

  std::vector<unsigned> key_tmp_vector ;

  const unsigned key_size = 2 + part_count ;

  unsigned * const key =
    ( key_size <= KEY_TMP_BUFFER_SIZE )
    ? key_tmp_buffer
    : ( key_tmp_vector.resize( key_size ) , & key_tmp_vector[0] );

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , family_count }
  // Thus family_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key.

  key[ key[0] = part_count + 1 ] = max ;

  {
    unsigned * const k = key + 1 ;
    for ( unsigned i = 0 ; i < part_count ; ++i ) { k[i] = part_ord[i] ; }
  }

  //----------------------------------
  // Bucket family has all of the same parts.
  // Look for the last bucket in this family:

  const std::vector<Bucket*>::iterator ik = lower_bound( bucket_set , key );

  //----------------------------------
  // If a member of the bucket family has space it is the last one
  // since buckets are kept packed.
  const bool bucket_family_exists =
    ik != bucket_set.begin() && bucket_part_equal( ik[-1]->key() , key );

  Bucket * const last_bucket = bucket_family_exists ? ik[-1] : NULL ;

  Bucket          * bucket    = NULL ;
  impl::BucketImpl::DataMap * field_map = NULL ;

  if ( last_bucket == NULL ) { // First bucket in this family
    key[ key[0] ] = 0 ; // Set the key's family count to zero
  }
  else { // Last bucket present, can it hold one more entity?

    if ( 0 == last_bucket->size() ) {
      throw std::logic_error( std::string(method) );
    }

    field_map = last_bucket->m_bucketImpl.m_field_map ;

    const unsigned last_count = last_bucket->key()[ key[0] ];

    const unsigned cap = last_bucket->capacity();

    if ( last_bucket->size() < cap ) {
      bucket = last_bucket ;
    }
    else if ( last_count < max ) {
      key[ key[0] ] = 1 + last_count ; // Increment the key's family count.
    }
    else {
      // ERROR insane number of buckets!
      std::string msg ;
      msg.append( method );
      msg.append( " FAILED due to insanely large number of buckets" );
      throw std::logic_error( msg );
    }
  }

  //----------------------------------
  // Family's field map does not exist, create it:

  if ( NULL == field_map ) {

    field_map = reinterpret_cast<impl::BucketImpl::DataMap*>(
                local_malloc( sizeof(impl::BucketImpl::DataMap) * ( num_fields + 1 )));

    // Start field data memory after the array of member entity pointers:
    unsigned value_offset = align( sizeof(Entity*) * bucket_capacity );

    for ( unsigned i = 0 ; i < num_fields ; ++i ) {
      const FieldBase  & field = * field_set[i] ;

      unsigned value_size = 0 ;

      const FieldBase::Restriction & dim =
        dimension( field, arg_entity_rank, part_count, part_ord, method);

      if ( dim.stride[0] ) { // Exists

        const unsigned type_stride = field.data_traits().stride_of ;
        const unsigned field_rank  = field.rank();

        value_size = type_stride *
          ( field_rank ? dim.stride[ field_rank - 1 ] : 1 );
      }

      field_map[i].m_base = value_offset ;
      field_map[i].m_size = value_size ;
      field_map[i].m_stride = dim.stride ;

      value_offset += align( value_size * bucket_capacity );
    }
    field_map[ num_fields ].m_base  = value_offset ;
    field_map[ num_fields ].m_size = 0 ;
    field_map[ num_fields ].m_stride = NULL ;
  }

  //----------------------------------

  if ( NULL == bucket ) {

    // Required bucket does not exist, must allocate and insert
    //
    // Allocation size:
    //   sizeof(Bucket) +
    //   key_size * sizeof(unsigned) +
    //   sizeof(Entity*) * capacity() +
    //   sum[number_of_fields]( fieldsize * capacity )
    //
    // The field_map[ num_fields ].m_base spans
    //   sizeof(Entity*) * capacity() +
    //   sum[number_of_fields]( fieldsize * capacity )

    const unsigned alloc_size = align( sizeof(Bucket) ) +
                                align( sizeof(unsigned) * key_size ) +
                                field_map[ num_fields ].m_base ;

    // All fields checked and sized, Ready to allocate

    void * const alloc_ptr = local_malloc( alloc_size );

    unsigned char * ptr = reinterpret_cast<unsigned char *>( alloc_ptr );

    ptr += align( sizeof( Bucket ) );

    unsigned * const new_key = reinterpret_cast<unsigned *>( ptr );

    ptr += align( sizeof(unsigned) * key_size );

    Entity ** const entity_array = reinterpret_cast<Entity**>( ptr );

    for ( unsigned i = 0 ; i < key_size ; ++i ) { new_key[i] = key[i] ; }

    bucket = new( alloc_ptr ) Bucket( mesh, arg_entity_rank , new_key,
                                      alloc_size, bucket_capacity ,
                                      field_map , entity_array );

    Bucket * first_bucket = last_bucket ? last_bucket->m_bucketImpl.m_bucket : bucket ;

    bucket->m_bucketImpl.m_bucket = first_bucket ; // Family members point to first bucket

    first_bucket->m_bucketImpl.m_bucket = bucket ; // First bucket points to new last bucket

    bucket_set.insert( ik , bucket );
  }

  //----------------------------------

  return bucket ;
}

//----------------------------------------------------------------------
// Every bucket in the family points to the first bucket,
// except the first bucket which points to the last bucket.

Bucket * BucketImpl::last_bucket_in_family( Bucket * b )
{
  static const char method[] = "stk::mesh::BucketImpl::last_bucket_in_family" ;

  Bucket * const first = bucket_counter( b->key() ) ? b->m_bucketImpl.m_bucket : b ;
  Bucket * const last  = first->m_bucketImpl.m_bucket ;

  if ( NULL == last || 0 == last->size() ) {
    throw std::logic_error( std::string(method) );
  }

  return last ;
}

//----------------------------------------------------------------------

void BucketImpl::zero_fields( Bucket & k_dst , unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const p = reinterpret_cast<unsigned char*>(k_dst.m_bucketImpl.m_entities);
  const DataMap *       i = k_dst.m_bucketImpl.m_field_map ;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i ) {
    if ( i->m_size ) {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void BucketImpl::copy_fields( Bucket & k_dst , unsigned i_dst ,
                          Bucket & k_src , unsigned i_src )
{
  static const char method[] = "stk::mesh::impl::BucketImpl::copy_fields" ;

  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const s = reinterpret_cast<unsigned char*>(k_src.m_bucketImpl.m_entities);
  unsigned char * const d = reinterpret_cast<unsigned char*>(k_dst.m_bucketImpl.m_entities);
  const DataMap *       j = k_src.m_bucketImpl.m_field_map ;
  const DataMap *       i = k_dst.m_bucketImpl.m_field_map ;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i , ++j ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        if ( i->m_size == j->m_size ) {
          memory_copy( d + i->m_base + i->m_size * i_dst ,
                       s + j->m_base + j->m_size * i_src , i->m_size );
        }
        else {
          std::ostringstream msg ;
          msg << method ;
          msg << " FAILED WITH INCOMPATIBLE FIELD SIZES" ;
          throw std::runtime_error( msg.str() );
        }
      }
      else {
        memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
      }
    }
  }
}


//----------------------------------------------------------------------

void BucketImpl::update_state()
{
  if ( 0 == bucket_counter( m_key ) ) {

    const MetaData & S = m_mesh.mesh_meta_data();
    const std::vector<FieldBase*> & field_set = S.get_fields();

    for ( unsigned i = 0 ; i < field_set.size() ; ) {

      DataMap * const tmp = m_field_map + i ;
      const FieldBase & field = * field_set[i] ;
      const unsigned num_state = field.number_of_states();
      i += num_state ;

      if ( 1 < num_state && tmp->m_size ) {
        unsigned offset[ MaximumFieldStates ] ;

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          offset[j] = tmp[j].m_base ;
        }

        for ( unsigned j = 0 ; j < num_state ; ++j ) {
          const unsigned j_new = ( j + num_state - 1 ) % num_state ;
          tmp[j_new].m_base = offset[j] ;
        }
      }
    }
  }
}



} // namespace impl 
} // namespace mesh
} // namespace stk



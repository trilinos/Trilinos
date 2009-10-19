#include <stdlib.h>
#include <memory.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace mesh {

namespace {

void memory_copy( unsigned char * dst , unsigned char * src , unsigned n )
{ memcpy( dst , src , n ); }

void memory_zero( unsigned char * dst , unsigned n )
{ memset( dst , 0 , n ); }

}

//----------------------------------------------------------------------
// BucketKey key = ( part-count , { part-ordinals } , counter )
//  key[ key[0] ] == counter

namespace {

inline
unsigned bucket_counter( const unsigned * const key )
{ return key[ *key ]; }

// The part count and parts are equal
bool bucket_part_equal( const unsigned * lhs , const unsigned * rhs )
{
  bool result = true ;
  {
    const unsigned * const end_lhs = lhs + *lhs ;
    while ( result && end_lhs != lhs ) {
      result = *lhs == *rhs ;
      ++lhs ; ++rhs ;
    }
  }
  return result ;
}

inline
bool bucket_key_less( const unsigned * lhs , const unsigned * rhs )
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

struct BucketLess {
  bool operator()( const Bucket * lhs_bucket , const unsigned * rhs ) const ;
  bool operator()( const unsigned * lhs , const Bucket * rhs_bucket ) const ;
};

// The part count and part ordinals are less
bool BucketLess::operator()( const Bucket * lhs_bucket ,
                             const unsigned * rhs ) const
{ return bucket_key_less( lhs_bucket->key() , rhs ); }

bool BucketLess::operator()( const unsigned * lhs ,
                             const Bucket * rhs_bucket ) const
{ return bucket_key_less( lhs , rhs_bucket->key() ); }

std::vector<Bucket*>::iterator
lower_bound( std::vector<Bucket*> & v , const unsigned * key )
{ return std::lower_bound( v.begin() , v.end() , key , BucketLess() ); }

}

//----------------------------------------------------------------------

bool Bucket::member( const Part & part ) const
{
  const unsigned * const i_beg = m_key + 1 ;
  const unsigned * const i_end = m_key + m_key[0] ;

  const unsigned ord = part.mesh_meta_data_ordinal();
  const unsigned * const i = std::lower_bound( i_beg , i_end , ord );

  return i_end != i && ord == *i ;
}

bool Bucket::member_all( const std::vector<Part*> & parts ) const
{
  const unsigned * const i_beg = m_key + 1 ;
  const unsigned * const i_end = m_key + m_key[0] ;

  const std::vector<Part*>::const_iterator ip_end = parts.end();
        std::vector<Part*>::const_iterator ip     = parts.begin() ; 

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_all = i_end != i && ord == *i ;
  }
  return result_all ;
}

bool Bucket::member_any( const std::vector<Part*> & parts ) const
{
  const unsigned * const i_beg = m_key + 1 ;
  const unsigned * const i_end = m_key + m_key[0] ;

  const std::vector<Part*>::const_iterator ip_end = parts.end();
        std::vector<Part*>::const_iterator ip     = parts.begin() ; 

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

//----------------------------------------------------------------------

bool has_superset( const Bucket & bucket , const Part & p )
{
  const unsigned ordinal = p.mesh_meta_data_ordinal();

  std::pair<const unsigned *, const unsigned *> 
    part_ord = bucket.superset_part_ordinals();

  part_ord.first =
    std::lower_bound( part_ord.first , part_ord.second , ordinal );

  return part_ord.first < part_ord.second && ordinal == *part_ord.first ;
}

bool has_superset( const Bucket & bucket , const PartVector & ps )
{
  const std::pair<const unsigned *, const unsigned *> 
    part_ord = bucket.superset_part_ordinals();

  bool result = ! ps.empty();

  for ( PartVector::const_iterator
        i = ps.begin() ; result && i != ps.end() ; ++i ) {

    const unsigned ordinal = (*i)->mesh_meta_data_ordinal();

    const unsigned * iter =
      std::lower_bound( part_ord.first , part_ord.second , ordinal );

    result = iter < part_ord.second && ordinal == *iter ;
  }
  return result ;
}

void Bucket::supersets( PartVector & ps ) const
{
  const MetaData & mesh_meta_data = m_mesh.mesh_meta_data();

  std::pair<const unsigned *, const unsigned *> 
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = & mesh_meta_data.get_part( * part_ord.first );
  }
}

//----------------------------------------------------------------------

bool field_data_valid( const FieldBase & f ,
                       const Bucket & k ,
                       unsigned ord ,
                       const char * required_by )
{
  const MetaData * const k_mesh_meta_data = & k.mesh().mesh_meta_data();
  const MetaData * const f_mesh_meta_data = & f.mesh_meta_data();
  const bool ok_mesh_meta_data  = k_mesh_meta_data == f_mesh_meta_data ;
  const bool ok_ord     = ord < k.size() ;
  const bool exists     = ok_mesh_meta_data && ok_ord &&
                          NULL != field_data( f , k.begin() );

  if ( required_by && ! exists ) {
    std::ostringstream msg ;
    msg << "stk::mesh::field_data_valid( " ;
    msg << f ;
    msg << " , " ;
    msg << k ;
    msg << " , " ;
    msg << ord ;
    msg << " , " ;
    msg << required_by ;
    msg << " ) FAILED with " ;
    if ( ! ok_mesh_meta_data ) {
      msg << " different MetaData" ;
    }
    else if ( ! ok_ord ) {
      msg << " Ordinal " ;
      msg << ord ;
      msg << " >= " ;
      msg << " size " ;
      msg << k.size();
    }
    else {
      msg << " no data" ;
    }
    throw std::runtime_error( msg.str() );
  }

  return exists ;
}

void throw_field_data_array( const FieldBase & f , unsigned R )
{
  std::ostringstream msg ;
  msg << "stk::mesh::throw_field_data_array( Field["
      << f.name() << "].rank() = " << f.rank()
      << " , truncation_rank = " << R
      << " ) FAILED, bad array truncation" ;
  throw std::runtime_error( msg.str() );
}

//----------------------------------------------------------------------

Bucket::Bucket( BulkData        & arg_mesh ,
                unsigned          arg_entity_type ,
                const unsigned  * arg_key ,
                size_t            arg_alloc_size ,
                size_t            arg_capacity ,
                Bucket::DataMap * arg_field_map ,
                Entity         ** arg_entity_array )
: m_mesh( arg_mesh ) ,
  m_entity_type( arg_entity_type ),
  m_key( arg_key ),
  m_alloc_size( arg_alloc_size ),
  m_capacity( arg_capacity ),
  m_size( 0 ),
  m_bucket(),
  m_field_map( arg_field_map ),
  m_entities( arg_entity_array )
{}


//----------------------------------------------------------------------

void Bucket::zero_fields( Bucket & k_dst , unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const p = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  const DataMap *       i = k_dst.m_field_map ;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i ) {
    if ( i->m_size ) {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void Bucket::copy_fields( Bucket & k_dst , unsigned i_dst ,
                          Bucket & k_src , unsigned i_src )
{
  static const char method[] = "stk::mesh::Bucket::copy_fields" ;

  const std::vector<FieldBase*> & field_set =
    k_dst.mesh().mesh_meta_data().get_fields();

  unsigned char * const s = reinterpret_cast<unsigned char*>(k_src.m_entities);
  unsigned char * const d = reinterpret_cast<unsigned char*>(k_dst.m_entities);
  const DataMap *       j = k_src.m_field_map ;
  const DataMap *       i = k_dst.m_field_map ;
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

}

//----------------------------------------------------------------------

void Bucket::update_state()
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

//----------------------------------------------------------------------

namespace {

void * local_malloc( size_t n )
{
  void * const ptr = malloc( n );

  if ( NULL == ptr ) {
    std::ostringstream msg ;
    msg << "stk::mesh::Bucket::declare_bucket FAILED malloc( " << n << " )" ;
    throw std::runtime_error( msg.str() );
  }

  return ptr ;
}

}

Bucket::~Bucket()
{}

void Bucket::destroy_bucket( Bucket * k )
{
  if ( 0 == bucket_counter( k->m_key ) ) {
    free( k->m_field_map );
  }

  k->~Bucket();

  free( k );
}

//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.

Bucket *
Bucket::declare_bucket( BulkData & mesh ,
                        const unsigned arg_entity_type ,
                        const unsigned part_count ,
                        const unsigned part_ord[] ,
                        const unsigned bucket_capacity ,
                        const std::vector< FieldBase * > & field_set ,
                              std::vector<Bucket*>       & bucket_set )
{
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  static const char method[] = "stk::mesh::Bucket::declare_bucket" ;

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

  const bool bucket_family_exists =
    ik != bucket_set.begin() && bucket_part_equal( ik[-1]->m_key , key );

  Bucket * const last_bucket = bucket_family_exists ? ik[-1] : NULL ;

  Bucket          * bucket    = NULL ;
  Bucket::DataMap * field_map = NULL ;

  if ( last_bucket == NULL ) { // First bucket in this family
    key[ key[0] ] = 0 ; // Set the key's family count to zero
  }
  else { // Last bucket present, can it hold one more entity?

    field_map = last_bucket->m_field_map ;

    const unsigned last_count = last_bucket->m_key[ key[0] ];

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

    field_map = reinterpret_cast<Bucket::DataMap*>(
                local_malloc( sizeof(Bucket::DataMap) * ( num_fields + 1 )));

    // Start field data memory after the array of member entity pointers:
    unsigned value_offset = align( sizeof(Entity*) * bucket_capacity );

    for ( unsigned i = 0 ; i < num_fields ; ++i ) {
      const FieldBase  & field = * field_set[i] ;

      unsigned value_size = 0 ;

      const FieldBase::Restriction & dim =
        dimension( field, arg_entity_type, part_count, part_ord, method);

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

    bucket = new( alloc_ptr ) Bucket( mesh, arg_entity_type , new_key,
                                      alloc_size, bucket_capacity ,
                                      field_map , entity_array );

    Bucket * first_bucket = last_bucket ? last_bucket->m_bucket : bucket ;

    bucket->m_bucket = first_bucket ; // Family members point to first bucket

    first_bucket->m_bucket = bucket ; // First bucket points to new last bucket

    bucket_set.insert( ik , bucket );
  }

  //----------------------------------

  return bucket ;
}

//----------------------------------------------------------------------
// The input part ordinals are complete and contain all supersets.

Bucket *
Bucket::declare_nil_bucket( BulkData & mesh , unsigned field_count )
{
  //----------------------------------
  // Field map gives NULL for all field data.

  Bucket::DataMap * field_map =
    reinterpret_cast<Bucket::DataMap*>(
      local_malloc( sizeof(Bucket::DataMap) * ( field_count + 1 )));

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

  const unsigned bad_entity_type = ~0u ;

  Bucket * bucket =
    new( alloc_ptr ) Bucket( mesh , bad_entity_type , new_key ,
                             alloc_size , 0 , field_map , NULL );

  bucket->m_bucket = bucket ;

  //----------------------------------

  return bucket ;
}

//----------------------------------------------------------------------

void BulkData::remove_entity( Bucket * k , unsigned i )
{
  Bucket * const first = bucket_counter( k->m_key ) ? k->m_bucket : k ;
  Bucket * const last  = first->m_bucket ;

  // Only move if not the last entity being removed

  if ( last != k || k->m_size != i + 1 ) {

    // Not the same bucket or not the last entity

    // Copy last entity in last to ik slot i

    Entity * const entity = last->m_entities[ last->m_size - 1 ];

    Bucket::copy_fields( *k , i , *last , last->m_size - 1 );

    k->m_entities[i]     = entity ;
    entity->m_bucket     = k ;
    entity->m_bucket_ord = i ;

    // Entity field data has relocated

    internal_propagate_relocation( *entity );
  }

  --( last->m_size );

  if ( last->m_size != 0 ) {
    last->m_entities[ last->m_size ] = NULL ;
  }
  else {

    // The current 'last' bucket is to be deleted.
    // The previous 'last' bucket becomes the
    // new 'last' bucket in the family:

    std::vector<Bucket*> & bucket_set = m_buckets[ last->entity_type() ];

    std::vector<Bucket*>::iterator ik = lower_bound(bucket_set, last->m_key);

    if ( ik == bucket_set.end() || last != *ik ) {
      throw std::runtime_error(
        std::string("stk::mesh::BulkData::remove_entity INTERNAL FAILURE") );
    }

    ik = bucket_set.erase( ik );

    if ( first != last ) { first->m_bucket = *--ik ; }

    Bucket::destroy_bucket( last );
  }
}

//----------------------------------------------------------------------

namespace {

struct LessEntityPointer {
  bool operator()( const Entity * const lhs , const Entity * const rhs ) const
    { return lhs->key() < rhs->key() ; }
};

}

void BulkData::internal_sort_bucket_entities()
{
  for ( unsigned entity_type = 0 ;
                 entity_type < m_buckets.size() ; ++entity_type ) {

    std::vector<Bucket*> & buckets = m_buckets[ entity_type ];

    size_t bk = 0 ; // Offset to first bucket of the family
    size_t ek = 0 ; // Offset to end   bucket of the family

    for ( ; bk < buckets.size() ; bk = ek ) {
      Bucket * ik_vacant = buckets[bk]->m_bucket ; // Last bucket, need space
      unsigned ie_vacant = ik_vacant->size();

      if ( ik_vacant->capacity() <= ie_vacant ) {
        // Have to create a bucket just for the scratch space...
        const unsigned * const bucket_key = buckets[bk]->m_key ;
        const unsigned         part_count = bucket_key[0] - 1 ;
        const unsigned * const part_ord   = bucket_key + 1 ;

        ik_vacant = Bucket::declare_bucket( *this , entity_type ,
                                            part_count , part_ord ,
                                            m_bucket_capacity ,
                                            m_mesh_meta_data.get_fields() ,
                                            buckets );
        ie_vacant = 0 ;
      }

      ik_vacant->m_entities[ ie_vacant ] = NULL ;

      // Determine offset to the end bucket in this family:
      while ( ek < buckets.size() && ik_vacant != buckets[ek] ) { ++ek ; }
      ++ek ;

      unsigned count = 0 ;
      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        count += buckets[ik]->size();
      }

      std::vector<Entity*> entities( count );

      std::vector<Entity*>::iterator j = entities.begin();

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          *j = b.m_entities[i] ;
        }
      }

      std::sort( entities.begin() , entities.end() , LessEntityPointer() );

      j = entities.begin();

      bool change_this_family = false ;

      for ( size_t ik = bk ; ik != ek ; ++ik ) {
        Bucket & b = * buckets[ik];
        const unsigned n = b.size();
        for ( unsigned i = 0 ; i < n ; ++i , ++j ) {
          Entity * const current = b.m_entities[i] ;

          if ( current != *j ) {

            if ( current ) {
              // Move current entity to the vacant spot
              Bucket::copy_fields( *ik_vacant , ie_vacant , b, i );
              current->m_bucket     = ik_vacant ;
              current->m_bucket_ord = ie_vacant ;
              ik_vacant->m_entities[ ie_vacant ] = current ;
            }

            // Set the vacant spot to where the required entity is now.
            ik_vacant = (*j)->m_bucket ;
            ie_vacant = (*j)->m_bucket_ord ;
            ik_vacant->m_entities[ ie_vacant ] = NULL ;

            // Move required entity to the required spot
            Bucket::copy_fields( b, i, *ik_vacant , ie_vacant );
            (*j)->m_bucket     = & b ;
            (*j)->m_bucket_ord = i ;
            b.m_entities[i]    = *j ;

            change_this_family = true ;
          }

          // Once a change has occured then need to propagate the
          // relocation for the remainder of the family.
          // This allows the propagation to be performed once per
          // entity as opposed to both times the entity is moved.

          if ( change_this_family ) { internal_propagate_relocation( **j ); }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Bucket & k )
{
  const MetaData & mesh_meta_data = k.mesh().mesh_meta_data();
  const std::string & entity_name =
    mesh_meta_data.entity_type_names()[ k.entity_type() ];

  PartVector parts ; k.supersets( parts );

  s << "Bucket( " << entity_name << " : " ;
  for ( PartVector::iterator i = parts.begin() ; i != parts.end() ; ++i ) {
    s << (*i)->name() << " " ;
  }
  s << ")" ;

  return s ;
}


std::ostream &
print( std::ostream & os , const std::string & indent , const Bucket & bucket )
{
  const MetaData & mesh_meta_data = bucket.mesh().mesh_meta_data();
  const std::string & entity_name =
    mesh_meta_data.entity_type_names()[ bucket.entity_type() ];

  const std::pair<const unsigned *, const unsigned *>
    part_ids = bucket.superset_part_ordinals();

  os << "Bucket(" << std::endl << indent << "Part intersection {" ;

  for ( const unsigned * i = part_ids.first ; i < part_ids.second ; ++i ) {
    const Part & part = mesh_meta_data.get_part( *i );
    os << " " << part.name();
  }

  os << " }" << std::endl << indent << entity_name << " members {" ;

  for ( unsigned j = 0 ; j < bucket.size() ; ++j ) {
    const EntityId id = bucket[j].identifier();
    os << " " << id ;
  }
  os << " } )" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

void
BucketIterator::throw_error(const char * err) const
{
  static const char header[] = "stk::mesh::Bucket::iterator FAILED: " ;
  std::string msg( header );
  msg.append( err );
  throw std::runtime_error( msg );
}


} // namespace mesh
} // namespace stk


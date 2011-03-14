/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------
#include <sstream>
#include <cstdlib>
#include <cstring>
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

void memory_copy( unsigned char * dst , const unsigned char * src , unsigned n )
{ std::memcpy( dst , src , n ); }


void memory_zero( unsigned char * dst , unsigned n )
{ std::memset( dst , 0 , n ); }

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
BucketImpl::~BucketImpl()
{
  bool this_is_first_bucket_in_family = (bucket_counter() == 0);
  if (this_is_first_bucket_in_family) {
    try {
      std::free( m_field_map );
    } catch(...) {}
  }
}

//----------------------------------------------------------------------

void BucketImpl::update_state()
{
  bool this_is_first_bucket_in_family = ( bucket_counter() == 0 );
  if (this_is_first_bucket_in_family) {

    const MetaData & S = MetaData::get(m_mesh);
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
// Every bucket in the family points to the first bucket,
// except the first bucket which points to the last bucket.

Bucket * BucketImpl::last_bucket_in_family()
{
  Bucket * last = last_bucket_in_family_impl();

  ThrowRequireMsg( NULL != last, "Last is NULL");
  ThrowRequireMsg( last->size() != 0, "Last bucket is empty");

  return last ;
}

Bucket * BucketImpl::last_bucket_in_family_impl()
{
  bool this_is_first_bucket_in_family = (bucket_counter() == 0);

  Bucket * last = NULL;

  if (this_is_first_bucket_in_family) {
    last = m_bucket;
  } else {
    last = m_bucket->m_bucketImpl.m_bucket;
  }

  return last;
}

//----------------------------------------------------------------------

Bucket * BucketImpl::first_bucket_in_family()
{
  return last_bucket_in_family_impl()->m_bucketImpl.m_bucket;
}

//----------------------------------------------------------------------

void BucketImpl::set_last_bucket_in_family( Bucket * last_bucket )
{
  Bucket * last = last_bucket_in_family_impl();
  Bucket * first = last->m_bucketImpl.m_bucket;
  first->m_bucketImpl.m_bucket = last_bucket;
}

//----------------------------------------------------------------------

void BucketImpl::set_first_bucket_in_family( Bucket * first_bucket )
{
  m_bucket = first_bucket;
}

//----------------------------------------------------------------------

BucketImpl::DataMap * BucketImpl::get_field_map()
{
  return m_field_map;
}

//----------------------------------------------------------------------

void BucketImpl::zero_fields( unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const p = reinterpret_cast<unsigned char*>(m_entities);
  const DataMap *       i = m_field_map;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i ) {
    if ( i->m_size ) {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void BucketImpl::replace_fields( unsigned i_dst , Bucket & k_src , unsigned i_src )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const s = reinterpret_cast<unsigned char*>(k_src.m_bucketImpl.m_entities);
  unsigned char * const d = reinterpret_cast<unsigned char*>(m_entities);
  const DataMap *       j = k_src.m_bucketImpl.m_field_map;
  const DataMap *       i = m_field_map;
  const DataMap * const e = i + field_set.size();

  for ( ; i != e ; ++i , ++j ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        ThrowErrorMsgIf( i->m_size != j->m_size,
            "Incompatible field sizes: " << i->m_size << " != " << j->m_size );

        memory_copy( d + i->m_base + i->m_size * i_dst ,
                     s + j->m_base + j->m_size * i_src , i->m_size );
      }
      else {
        memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
      }
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk



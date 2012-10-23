/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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
#include <stk_mesh/base/FindRestriction.hpp>

//----------------------------------------------------------------------
namespace {
inline unsigned align( size_t nb )
{
  enum { BYTE_ALIGN = 16 };
  const unsigned gap = nb % BYTE_ALIGN ;
  if ( gap ) { nb += BYTE_ALIGN - gap ; }
  return nb ;
}
}//namespace anonymous

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

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

// The part count and part ordinals are less
bool BucketLess::operator()( const Bucket * lhs_bucket ,
                             const unsigned * rhs ) const
{ return bucket_key_less( lhs_bucket->key() , rhs ); }

bool BucketLess::operator()( const unsigned * lhs ,
                             const Bucket * rhs_bucket ) const
{ return bucket_key_less( lhs , rhs_bucket->key() ); }

//----------------------------------------------------------------------

Bucket::Bucket( BulkData & arg_mesh ,
                EntityRank arg_entity_rank,
                const std::vector<unsigned> & arg_key,
                size_t arg_capacity
        )
  : m_mesh(arg_mesh)
  , m_entity_rank(arg_entity_rank)
  , m_key(arg_key)
  , m_capacity(arg_capacity)
  , m_size(0)
  , m_bucket(NULL)
  , m_field_map( m_mesh.mesh_meta_data().get_fields().size()+1)
  , m_entities(arg_capacity)
  , m_field_data(NULL)
  , m_field_data_end(NULL)
  , m_bucket_family(NULL)
{
  //calculate the size of the field_data

  const std::vector< FieldBase * > & field_set =
    arg_mesh.mesh_meta_data().get_fields();

  const size_t num_fields = field_set.size();

  size_t field_data_size = 0;

  if (arg_capacity != 0) {
    for ( size_t i = 0; i<num_fields; ++i) {
      const FieldBase  & field = * field_set[i] ;
      unsigned num_bytes_per_entity = 0 ;

      const FieldBase::Restriction & restriction =
        find_restriction( field, arg_entity_rank, &m_key[1], &m_key[1]+(m_key[0]-1), PartOrdLess());

      if ( restriction.dimension() > 0 ) { // Exists

        const unsigned type_stride = field.data_traits().stride_of ;
        const unsigned field_rank  = field.rank();

        num_bytes_per_entity = type_stride *
          ( field_rank ? restriction.stride( field_rank - 1 ) : 1 );
      }
      m_field_map[i].m_base = field_data_size ;
      m_field_map[i].m_size = num_bytes_per_entity ;
      m_field_map[i].m_stride = &restriction.stride(0);

      field_data_size += align( num_bytes_per_entity * m_capacity );
    }

    m_field_map[ num_fields ].m_base  = field_data_size ;
    m_field_map[ num_fields ].m_size = 0 ;
    m_field_map[ num_fields ].m_stride = NULL ;
  }
  else { //nil bucket

    FieldBase::Restriction::size_type empty_stride[ MaximumFieldDimension ];
    Copy<MaximumFieldDimension>( empty_stride , FieldBase::Restriction::size_type(0) );

    for ( size_t i = 0; i<num_fields; ++i) {
      m_field_map[i].m_base = 0 ;
      m_field_map[i].m_size = 0 ;
      m_field_map[i].m_stride = empty_stride;
    }
    m_field_map[ num_fields ].m_base   = 0 ;
    m_field_map[ num_fields ].m_size   = 0 ;
    m_field_map[ num_fields ].m_stride = NULL ;
  }

  //allocate space for the fields
  m_field_data = field_data_size > 0 ? new unsigned char[field_data_size] : NULL;

  //
  //[TODO] ALAN, TODD: to investigate if memory_zero is necessary to fix valgrind
  //issues in the following regression test:
  //adagio_rtest/presto/super_elem_rigid_body/super_elem_rigid_body.test|np1_explicit_reverseMPC
  //memory_zero(m_field_data, field_data_size);
  //ABW UPDATE: found and fixed bug in strumento/src/element/SuperElementHandler.C
  //which was reading past end of field, so this memory_zero is no longer
  //necessary. 8/9/2012
  //std::memset( m_field_data , 0xfff , field_data_size );

  m_field_data_end = m_field_data + field_data_size;
}




bool Bucket::member( const Part & part ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const unsigned ord = part.mesh_meta_data_ordinal();
  const unsigned * const i = std::lower_bound( i_beg , i_end , ord );

  return i_end != i && ord == *i ;
}

bool Bucket::member_all( const PartVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_all = i_end != i && ord == *i ;
  }
  return result_all ;
}

bool Bucket::member_any( const PartVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

bool Bucket::member_any( const OrdinalVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const OrdinalVector::const_iterator ip_end = parts.end();
        OrdinalVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = *ip;
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

//----------------------------------------------------------------------

bool has_superset( const Bucket & bucket, const unsigned & ordinal )
{
  std::pair<const unsigned *, const unsigned *>
    part_ord = bucket.superset_part_ordinals();

  part_ord.first =
    std::lower_bound( part_ord.first , part_ord.second , ordinal );

  return part_ord.first < part_ord.second && ordinal == *part_ord.first ;
}

bool has_superset( const Bucket & bucket , const Part & p )
{
  const unsigned ordinal = p.mesh_meta_data_ordinal();
  return has_superset(bucket,ordinal);
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
  const MetaData & mesh_meta_data = MetaData::get( *this );

  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = & mesh_meta_data.get_part( * part_ord.first );
  }
}

void Bucket::supersets( OrdinalVector & ps ) const
{
  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = *part_ord.first;
  }
}

//----------------------------------------------------------------------

bool field_data_valid( const FieldBase & f ,
                       const Bucket & k ,
                       unsigned ord ,
                       const char * required_by )
{
  const MetaData * const k_mesh_meta_data = & MetaData::get(k);
  const MetaData * const f_mesh_meta_data = & MetaData::get(f);
  const bool ok_mesh_meta_data  = k_mesh_meta_data == f_mesh_meta_data ;
  const bool ok_ord     = ord < k.size() ;
  const bool exists     = ok_mesh_meta_data && ok_ord &&
                          NULL != field_data( f , k.begin() );

  if ( required_by && ! exists ) {
    std::ostringstream msg_begin ;
    msg_begin << "For args: " ;
    msg_begin << f << " , " ;
    msg_begin << k << " , " ;
    msg_begin << ord << " , " ;
    msg_begin << required_by ;
    msg_begin << "; operation FAILED with " ;
    ThrowErrorMsgIf( ! ok_mesh_meta_data,
                     msg_begin.str() << " different MetaData");
    ThrowErrorMsgIf( ! ok_ord, msg_begin.str() <<
                     " Ordinal " <<  ord << " >= " << " size " << k.size());
    ThrowErrorMsg( msg_begin.str() << " no data");
  }

  return exists ;
}

//----------------------------------------------------------------------
bool Bucket::assert_correct() const {
  // test equivalent() method
  const Bucket* bucket = this;
  const Bucket * first = first_bucket_in_family();
  if (!first || ! bucket->equivalent(*first) || ! first->equivalent(*bucket) )
    return false;

  // other tests...

  return true;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Bucket & k )
{
  const MetaData & mesh_meta_data = MetaData::get(k);
  const std::string & entity_rank_name =
    k.entity_rank() == InvalidEntityRank ? "Nil" :
    mesh_meta_data.entity_rank_names()[ k.entity_rank() ];

  PartVector parts ; k.supersets( parts );

  s << "Bucket( " << entity_rank_name << " : " ;
  for ( PartVector::iterator i = parts.begin() ; i != parts.end() ; ++i ) {
    s << (*i)->name() << " " ;
  }
  s << ")" ;

  return s ;
}


std::ostream &
print( std::ostream & os , const std::string & indent , const Bucket & bucket )
{
  const MetaData & mesh_meta_data = MetaData::get(bucket);
  const std::string & entity_rank_name =
    bucket.entity_rank() == InvalidEntityRank ? "Nil" :
    mesh_meta_data.entity_rank_names()[ bucket.entity_rank() ];

  const std::pair<const unsigned *, const unsigned *>
    part_ids = bucket.superset_part_ordinals();

  os << "Bucket(" << std::endl << indent << "Part intersection {" ;

  for ( const unsigned * i = part_ids.first ; i < part_ids.second ; ++i ) {
    const Part & part = mesh_meta_data.get_part( *i );
    os << " " << part.name();
  }

  os << " }" << std::endl << indent << entity_rank_name << " members {" ;

  for ( unsigned j = 0 ; j < bucket.size() ; ++j ) {
    const EntityId id = bucket[j].identifier();
    os << " " << id ;
  }
  os << " } )" << std::endl ;

  return os ;
}


namespace {

void memory_copy( unsigned char * dst , const unsigned char * src , unsigned n )
{ std::memcpy( dst , src , n ); }


void memory_zero( unsigned char * dst , unsigned n )
{ std::memset( dst , 0 , n ); }

} // namespace

//----------------------------------------------------------------------

void Bucket::update_state()
{
  const MetaData & meta = MetaData::get(m_mesh);
  const std::vector<FieldBase*> & field_set = meta.get_fields();

  for ( unsigned i = 0 ; i < field_set.size() ; ) {

    DataMap * const tmp = &m_field_map[0] + i ;
    const FieldBase & field = * field_set[i] ;
    const unsigned num_state = field.number_of_states();
    i += num_state ;

    if ( 1 < num_state && tmp->m_size ) {
      unsigned offset[ MaximumFieldStates ] ;

      offset[0] = tmp[num_state-1].m_base;
      for ( unsigned j = 1 ; j < num_state ; ++j ) {
        offset[j] = tmp[j-1].m_base ;
      }

      for ( unsigned j = 0 ; j < num_state ; ++j ) {
        tmp[j].m_base = offset[j] ;
      }
    }
  }
}

//----------------------------------------------------------------------
// Every bucket in the family points to the first bucket,
// except the first bucket which points to the last bucket.

Bucket * Bucket::last_bucket_in_family() const
{
  Bucket * last = last_bucket_in_family_impl();

  ThrowRequireMsg( NULL != last, "Last is NULL");
  ThrowRequireMsg( last->size() != 0, "Last bucket is empty");

  return last ;
}

Bucket * Bucket::last_bucket_in_family_impl() const
{
  bool this_is_first_bucket_in_family = (bucket_counter() == 0);

  Bucket * last = NULL;

  if (this_is_first_bucket_in_family) {
    last = m_bucket;
  } else {
    last = m_bucket->m_bucket;
  }

  return last;
}

//----------------------------------------------------------------------

Bucket * Bucket::first_bucket_in_family() const
{
  return last_bucket_in_family_impl()->m_bucket;
}

//----------------------------------------------------------------------

void Bucket::set_last_bucket_in_family( Bucket * last_bucket )
{
  Bucket * last = last_bucket_in_family_impl();
  Bucket * first = last->m_bucket;
  first->m_bucket = last_bucket;
}

//----------------------------------------------------------------------

void Bucket::set_first_bucket_in_family( Bucket * first_bucket )
{
  m_bucket = first_bucket;
}

//----------------------------------------------------------------------

Bucket::DataMap * Bucket::get_field_map()
{
  return &m_field_map[0];
}

//----------------------------------------------------------------------

void Bucket::initialize_fields( unsigned i_dst )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const p = m_field_data;
  const DataMap *       i = &m_field_map[0];
  const DataMap * const e = i + field_set.size();

  for (std::vector<FieldBase*>::const_iterator field_iter=field_set.begin() ;
       i != e ; ++i, ++field_iter ) {

    if (i->m_size == 0) {
      continue;
    }

    const unsigned char* init_val = reinterpret_cast<const unsigned char*>((*field_iter)->get_initial_value());
    if (init_val != NULL) {
      memory_copy( p + i->m_base + i->m_size * i_dst , init_val, i->m_size );
    }
    else {
      memory_zero( p + i->m_base + i->m_size * i_dst , i->m_size );
    }
  }
}

void Bucket::replace_fields( unsigned i_dst , Bucket & k_src , unsigned i_src )
{
  const std::vector<FieldBase*> & field_set =
    MetaData::get(m_mesh).get_fields();

  unsigned char * const s = k_src.m_field_data;
  unsigned char * const d = m_field_data;
  const DataMap *       j = &(k_src.m_field_map[0]);
  const DataMap *       i = &m_field_map[0];
  const DataMap * const e = i + field_set.size();

  for (std::vector<FieldBase*>::const_iterator field_iter=field_set.begin() ;
       i != e ; ++i , ++j, ++field_iter ) {

    if ( i->m_size ) {
      if ( j->m_size ) {
        ThrowAssertMsg( i->m_size == j->m_size,
            "Incompatible field sizes: " << i->m_size << " != " << j->m_size );

        memory_copy( d + i->m_base + i->m_size * i_dst ,
                     s + j->m_base + j->m_size * i_src , i->m_size );
      }
      else {
        const unsigned char* init_val = reinterpret_cast<const unsigned char*>((*field_iter)->get_initial_value());
        if (init_val != NULL) {
          memory_copy( d + i->m_base + i->m_size * i_dst ,
                       init_val, i->m_size );
        }
        else {
          memory_zero( d + i->m_base + i->m_size * i_dst , i->m_size );
        }
      }
    }
  }
}



} // namespace mesh
} // namespace stk


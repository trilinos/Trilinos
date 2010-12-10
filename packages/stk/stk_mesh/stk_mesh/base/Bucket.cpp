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

bool Bucket::member( const Part & part ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const unsigned ord = part.mesh_meta_data_ordinal();
  const unsigned * const i = std::lower_bound( i_beg , i_end , ord );

  return i_end != i && ord == *i ;
}

bool Bucket::member_all( const std::vector<Part*> & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

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
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

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
  const MetaData & mesh_meta_data = mesh().mesh_meta_data();

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
                EntityRank        arg_entity_rank ,
                const unsigned  * arg_key ,
                size_t            arg_alloc_size ,
                size_t            arg_capacity ,
                impl::BucketImpl::DataMap * arg_field_map ,
                Entity         ** arg_entity_array )
: m_bucketImpl( arg_mesh, arg_entity_rank, arg_key, arg_alloc_size, arg_capacity, arg_field_map, arg_entity_array )
{}


//----------------------------------------------------------------------

Bucket::~Bucket()
{}



//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Bucket & k )
{
  const MetaData & mesh_meta_data = k.mesh().mesh_meta_data();
  const std::string & entity_name =
    mesh_meta_data.entity_rank_names()[ k.entity_rank() ];

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
    mesh_meta_data.entity_rank_names()[ bucket.entity_rank() ];

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


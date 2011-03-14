/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>

#include <cstring>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/util/SimpleArrayOps.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------

namespace {

struct FieldRestrictionLess {
  bool operator()( const FieldRestriction & lhs ,
                   const FieldRestriction & rhs ) const
    { return lhs.key < rhs.key ; }

  bool operator()( const FieldRestriction & lhs ,
                   const EntityKey & rhs ) const
    { return lhs.key < rhs ; }
};

FieldRestrictionVector::const_iterator
  find( const FieldRestrictionVector & v , const EntityKey & key )
{
  FieldRestrictionVector::const_iterator
    i = std::lower_bound( v.begin() , v.end() , key , FieldRestrictionLess() );

  if ( i != v.end() && i->key != key ) { i = v.end(); }

  return i ;
}

}

//----------------------------------------------------------------------

FieldBaseImpl::FieldBaseImpl(
    MetaData                   * arg_mesh_meta_data ,
    unsigned                     arg_ordinal ,
    const std::string          & arg_name ,
    const DataTraits           & arg_traits ,
    unsigned                     arg_rank,
    const shards::ArrayDimTag  * const * arg_dim_tags,
    unsigned                     arg_number_of_states ,
    FieldState                   arg_this_state
    )
: m_name( arg_name ),
  m_attribute(),
  m_data_traits( arg_traits ),
  m_meta_data( arg_mesh_meta_data ),
  m_ordinal( arg_ordinal ),
  m_num_states( arg_number_of_states ),
  m_this_state( arg_this_state ),
  m_rank( arg_rank ),
  m_dim_map()
{
  FieldBase * const pzero = NULL ;
  const shards::ArrayDimTag * const dzero = NULL ;
  Copy<MaximumFieldStates>(    m_field_states , pzero );
  Copy<MaximumFieldDimension>( m_dim_tags ,     dzero );

  for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
    m_dim_tags[i] = arg_dim_tags[i];
  }
}

const FieldRestrictionVector & FieldBaseImpl::restrictions() const
{ return m_field_states[0]->m_impl.m_dim_map ; }

FieldRestrictionVector & FieldBaseImpl::restrictions()
{ return m_field_states[0]->m_impl.m_dim_map ; }


//----------------------------------------------------------------------

namespace {

void print_restriction( std::ostream & os ,
                        unsigned type ,
                        const Part & part ,
                        unsigned rank ,
                        const FieldRestriction::size_type * stride )
{
  os << "{ entity_rank(" << type << ") part(" << part.name() << ") : " ;
  os << stride[0] ;
  for ( unsigned i = 1 ; i < rank ; ++i ) {
    if ( ! stride[i] ) {
      os << " , 0 " ;
    }
    else if ( stride[i] % stride[i-1] ) {
      os << " , " << stride[i] << " / " << stride[i-1] ;
    }
    else {
      os << " , " << stride[i] / stride[i-1] ;
    }
  }
  os << " }" ;
}

std::string print_restriction( unsigned type ,
                               const Part & part ,
                               unsigned rank ,
                               const FieldRestriction::size_type * stride )
{
  std::ostringstream oss;
  print_restriction(oss, type, part, rank, stride);
  return oss.str();
}

}

// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void FieldBaseImpl::insert_restriction(
  const char     * arg_method ,
  EntityRank       arg_entity_rank ,
  const Part     & arg_part ,
  const unsigned * arg_stride )
{
  FieldRestriction tmp ;

  tmp.key = EntityKey( arg_entity_rank , arg_part.mesh_meta_data_ordinal() );

  {
    unsigned i = 0 ;
    if ( m_rank ) {
      for ( i = 0 ; i < m_rank ; ++i ) { tmp.stride[i] = arg_stride[i] ; }
    }
    else { // Scalar field is 0 == m_rank
      i = 1 ;
      tmp.stride[0] = 1 ;
    }
    // Remaining dimensions are 1, no change to stride
    for ( ; i < MaximumFieldDimension ; ++i ) {
      tmp.stride[i] = tmp.stride[i-1] ;
    }

    for ( i = 1 ; i < m_rank ; ++i ) {
      const bool bad_stride = 0 == tmp.stride[i] ||
                              0 != tmp.stride[i] % tmp.stride[i-1];
      ThrowErrorMsgIf( bad_stride,
          arg_method << " FAILED for " << *this <<
          " WITH BAD STRIDE " <<
          print_restriction( arg_entity_rank, arg_part, m_rank, tmp.stride));;
    }
  }

  {
    FieldRestrictionVector & rMap = restrictions();

    FieldRestrictionVector::iterator i = rMap.begin(), j = rMap.end();

    i = std::lower_bound(i,j,tmp,FieldRestrictionLess());

    if ( i == j || i->key != tmp.key ) {
      rMap.insert( i , tmp );
    }
    else {
      ThrowErrorMsgIf( Compare<MaximumFieldDimension>::not_equal(i->stride,tmp.stride),
          arg_method << " FAILED for " << *this << " " <<
          print_restriction( arg_entity_rank, arg_part, m_rank, i->stride ) <<
          " WITH INCOMPATIBLE REDECLARATION " <<
          print_restriction( arg_entity_rank, arg_part, m_rank, tmp.stride ));
    }
  }
}

void FieldBaseImpl::verify_and_clean_restrictions(
  const char       * arg_method ,
  const PartVector & arg_all_parts )
{
  const EntityKey invalid_key ;
  FieldRestrictionVector & rMap = restrictions();
  FieldRestrictionVector::iterator i , j ;

  for ( i = rMap.begin() ; i != rMap.end() ; ++i ) {
    if ( i->key != invalid_key ) {
      const unsigned typeI = entity_rank( i->key );
      const Part   & partI = * arg_all_parts[ entity_id( i->key ) ];
      bool  found_superset = false ;

      for ( j = i + 1 ; j != rMap.end() && ! found_superset ; ++j ) {
        if ( j->key != invalid_key ) {
          const unsigned typeJ = entity_rank( j->key );
          const Part   & partJ = * arg_all_parts[ entity_id( j->key ) ];

          if ( typeI == typeJ ) {
            const bool found_subset = contain( partI.subsets() , partJ );
            found_superset = ! found_subset &&
                             contain( partI.supersets() , partJ );

            if ( found_subset || found_superset ) {
              ThrowErrorMsgIf( Compare< MaximumFieldDimension >::not_equal( i->stride , j->stride ),
                  arg_method << "[" << *this << "] FAILED: " <<
                  print_restriction( typeI, partI, m_rank, i->stride ) <<
                  ( found_subset ? " INCOMPATIBLE SUBSET " : " INCOMPATIBLE SUPERSET ") <<
                  print_restriction( typeJ, partJ, m_rank, j->stride ));
            }

            if ( found_subset ) { j->key = invalid_key; }
          }
        }
        if ( found_superset ) { i->key = invalid_key; }
      }
    }
  }

  // Clean out redundant entries:
  // NOTE: test for 'i != rMap.end()' not needed since i is
  //       incremented no more than j and j is checked. Keeping check in
  //       silences klocwork and guards against future change...
  for ( j = i = rMap.begin() ; j != rMap.end() && i != rMap.end(); ++j ) {
    if ( j->key != invalid_key ) {
      if ( i->key == invalid_key ) {
        *i = *j ;
      }
      ++i ;
    }
  }

  rMap.erase( i , j );
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// This part or any superset of this part

const FieldRestriction &
FieldBaseImpl::restriction( unsigned entity_rank , const Part & part ) const
{
  static const FieldRestriction empty ;

  const FieldRestrictionVector & rMap = restrictions();
  const FieldRestrictionVector::const_iterator ie = rMap.end() ;
        FieldRestrictionVector::const_iterator i ;

  const PartVector::const_iterator ipe = part.supersets().end();
        PartVector::const_iterator ip  = part.supersets().begin() ;

  // Start with this part:
  EntityKey key = EntityKey( entity_rank , part.mesh_meta_data_ordinal() );

  while ( ie == ( i = find( rMap , key ) ) && ipe != ip ) {
    // Not found try another superset part:
    key = EntityKey( entity_rank , (*ip)->mesh_meta_data_ordinal() );
    ++ip ;
  }

  return ie == i ? empty : *i ;
}

unsigned FieldBaseImpl::max_size( unsigned entity_rank ) const
{
  unsigned max = 0 ;

  const FieldRestrictionVector & rMap = restrictions();
  const FieldRestrictionVector::const_iterator ie = rMap.end() ;
        FieldRestrictionVector::const_iterator i = rMap.begin();

  for ( ; i != ie ; ++i ) {
    if ( i->type() == entity_rank ) {
      const unsigned len = m_rank ? i->stride[ m_rank - 1 ] : 1 ;
      if ( max < len ) { max = len ; }
    }
  }

  return max ;
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const FieldBaseImpl & field )
{
  s << "FieldBaseImpl<" ;
  s << field.data_traits().name ;
  for ( unsigned i = 0 ; i < field.rank() ; ++i ) {
    s << "," << field.dimension_tags()[i]->name();
  }
  s << ">" ;

  s << "[ name = \"" ;
  s << field.name() ;
  s << "\" , #states = " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  const PartVector & all_parts = MetaData::get(field).get_parts();
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();
  s << field.name() ;
  s << " {" ;
  for ( FieldBase::RestrictionVector::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << std::endl << b << "  " ;
    print_restriction( s, entity_rank( i->key ),
                       * all_parts[ entity_id( i->key ) ],
	               field.rank(), i->stride);
    s << std::endl;
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------




} // namespace impl
} // namespace mesh
} // namespace stk

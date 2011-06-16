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
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------

namespace {

FieldRestrictionVector::const_iterator
  find( const FieldRestrictionVector & v , const FieldRestriction & restr )
{
  FieldRestrictionVector::const_iterator
    i = std::lower_bound( v.begin() , v.end() , restr );

  if ( i != v.end() && !(*i == restr) ) { i = v.end(); }

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
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::FieldBaseImpl", LOG_FIELD, m_ordinal);

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
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::insert_restriction", LOG_FIELD, m_ordinal);

  FieldRestriction tmp( arg_entity_rank , arg_part.mesh_meta_data_ordinal() );

  {
    unsigned i = 0 ;
    if ( m_rank ) {
      for ( i = 0 ; i < m_rank ; ++i ) { tmp.stride(i) = arg_stride[i] ; }
    }
    else { // Scalar field is 0 == m_rank
      i = 1 ;
      tmp.stride(0) = 1 ;
    }
    // Remaining dimensions are 1, no change to stride
    for ( ; i < MaximumFieldDimension ; ++i ) {
      tmp.stride(i) = tmp.stride(i-1) ;
    }

    for ( i = 1 ; i < m_rank ; ++i ) {
      const bool bad_stride = 0 == tmp.stride(i) ||
                              0 != tmp.stride(i) % tmp.stride(i-1);
      ThrowErrorMsgIf( bad_stride,
          arg_method << " FAILED for " << *this <<
          " WITH BAD STRIDE " <<
          print_restriction( tmp, arg_entity_rank, arg_part, m_rank ));;
    }
  }

  {
    FieldRestrictionVector & rMap = restrictions();

    FieldRestrictionVector::iterator i = rMap.begin(), j = rMap.end();

    i = std::lower_bound(i,j,tmp);

    if ( i == j || !(*i == tmp) ) {
      rMap.insert( i , tmp );
    }
    else {
      ThrowErrorMsgIf( i->not_equal_stride(tmp),
          arg_method << " FAILED for " << *this << " " <<
          print_restriction( *i, arg_entity_rank, arg_part, m_rank ) <<
          " WITH INCOMPATIBLE REDECLARATION " <<
          print_restriction( tmp, arg_entity_rank, arg_part, m_rank ));
    }
  }
}

void FieldBaseImpl::verify_and_clean_restrictions(
  const char       * arg_method ,
  const PartVector & arg_all_parts )
{
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::verify_and_clean_restrictions", LOG_FIELD, m_ordinal);

  const FieldRestriction invalid_restr ;
  FieldRestrictionVector & rMap = restrictions();

  // Search for redundant field restrictions. A restriction R is redundant if there exists
  // another restriction on a superset of R's part that is for entities of the same rank as R.
  for (FieldRestrictionVector::iterator i = rMap.begin() ; i != rMap.end() ; ++i ) {
    if ( *i != invalid_restr ) {
      const EntityRank rankI = i->entity_rank();
      const Part     & partI = * arg_all_parts[ i->part_ordinal() ];
      bool    found_superset = false ;

      for (FieldRestrictionVector::iterator j = i + 1 ; j != rMap.end() && ! found_superset ; ++j ) {
        if ( *j != invalid_restr ) {
          const EntityRank rankJ = j->entity_rank();
          const Part     & partJ = * arg_all_parts[ j->part_ordinal() ];

          if ( rankI == rankJ ) {
            const bool found_subset = contain( partI.subsets() , partJ );
            found_superset = ! found_subset &&
                             contain( partI.supersets() , partJ );

            if ( found_subset || found_superset ) {
              ThrowErrorMsgIf( i->not_equal_stride(*j),
                  arg_method << "[" << *this << "] FAILED: " <<
                  print_restriction( *i, rankI, partI, m_rank ) <<
                  ( found_subset ? " INCOMPATIBLE SUBSET " : " INCOMPATIBLE SUPERSET ") <<
                  print_restriction( *j, rankJ, partJ, m_rank ));
            }

            if ( found_subset ) { *j = invalid_restr; }
          }
        }
        if ( found_superset ) { *i = invalid_restr; }
      }
    }
  }

  // Clean out redundant entries:
  FieldRestrictionVector::iterator new_end = std::remove(rMap.begin(), rMap.end(), invalid_restr);
  rMap.erase(new_end, rMap.end());
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
  FieldRestriction restr( entity_rank , part.mesh_meta_data_ordinal() );

  while ( ie == ( i = find( rMap , restr ) ) && ipe != ip ) {
    // Not found try another superset part:
    restr = FieldRestriction( entity_rank , (*ip)->mesh_meta_data_ordinal() );
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
    if ( i->entity_rank() == entity_rank ) {
      const unsigned len = m_rank ? i->stride( m_rank - 1 ) : 1 ;
      if ( max < len ) { max = len ; }
    }
  }

  return max ;
}

void FieldBaseImpl::set_field_states( FieldBase ** field_states)
{
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::set_field_states", LOG_FIELD, m_ordinal);

  for (unsigned i = 0; i < m_num_states; ++i) {
    m_field_states[i] = field_states[i];
  }
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
    i->print( s, i->entity_rank(), * all_parts[ i->part_ordinal() ], field.rank() );
    s << std::endl;
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------




} // namespace impl
} // namespace mesh
} // namespace stk

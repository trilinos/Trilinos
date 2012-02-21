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
  m_field_rank( arg_rank ),
  m_dim_map(),
  m_selector_restrictions(),
  m_initial_value(NULL),
  m_initial_value_num_bytes(0)
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

//----------------------------------------------------------------------
FieldBaseImpl::~FieldBaseImpl()
{
  if (state() == StateNone) {
    void*& init_val = m_initial_value;

    delete [] reinterpret_cast<char*>(init_val);
    init_val = NULL;
  }
}

//----------------------------------------------------------------------
const FieldRestrictionVector & FieldBaseImpl::restrictions() const
{ return m_field_states[0]->m_impl.m_dim_map ; }

const FieldRestrictionVector & FieldBaseImpl::selector_restrictions() const
{ return m_field_states[0]->m_impl.m_selector_restrictions ; }

FieldRestrictionVector & FieldBaseImpl::restrictions()
{ return m_field_states[0]->m_impl.m_dim_map ; }

FieldRestrictionVector & FieldBaseImpl::selector_restrictions()
{ return m_field_states[0]->m_impl.m_selector_restrictions ; }


//----------------------------------------------------------------------

// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void FieldBaseImpl::insert_restriction(
  const char     * arg_method ,
  EntityRank       arg_entity_rank ,
  const Part     & arg_part ,
  const unsigned * arg_stride,
  const void*      arg_init_value )
{
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::insert_restriction", LOG_FIELD, m_ordinal);

  FieldRestriction tmp( arg_entity_rank , arg_part.mesh_meta_data_ordinal() );

  {
    unsigned i = 0 ;
    if ( m_field_rank ) {
      for ( i = 0 ; i < m_field_rank ; ++i ) { tmp.stride(i) = arg_stride[i] ; }
    }
    else { // Scalar field is 0 == m_field_rank
      i = 1 ;
      tmp.stride(0) = 1 ;
    }
    // Remaining dimensions are 1, no change to stride
    for ( ; i < MaximumFieldDimension ; ++i ) {
      tmp.stride(i) = tmp.stride(i-1) ;
    }

    for ( i = 1 ; i < m_field_rank ; ++i ) {
      const bool bad_stride = 0 == tmp.stride(i) ||
                              0 != tmp.stride(i) % tmp.stride(i-1);
      ThrowErrorMsgIf( bad_stride,
          arg_method << " FAILED for " << *this <<
          " WITH BAD STRIDE " <<
          print_restriction( tmp, arg_entity_rank, arg_part, m_field_rank ));;
    }
  }

  if (arg_init_value != NULL) {
    //insert_restriction can be called multiple times for the same field, giving
    //the field different lengths on different mesh-parts.
    //We will only store one initial-value array, we need to store the one with
    //maximum length for this field so that it can be used to initialize data
    //for all field-restrictions. For the parts on which the field is shorter,
    //a subset of the initial-value array will be used.
    //
    //We want to end up storing the longest arg_init_value array for this field.
    //
    //Thus, we call set_initial_value only if the current length is longer
    //than what's already been stored.

    //length in bytes is num-scalars X sizeof-scalar:

    size_t num_scalars = 1;
    //if rank > 0, then field is not a scalar field, so num-scalars is
    //obtained from the stride array:
    if (m_field_rank > 0) num_scalars = tmp.stride(m_field_rank-1);

    size_t sizeof_scalar = m_data_traits.size_of;
    size_t nbytes = sizeof_scalar * num_scalars;

    size_t old_nbytes = 0;
    if (get_initial_value() != NULL) {
      old_nbytes = get_initial_value_num_bytes();
    }

    if (nbytes > old_nbytes) {
      set_initial_value(arg_init_value, num_scalars, nbytes);
    }
  }

  {
    FieldRestrictionVector & restrs = restrictions();

    FieldRestrictionVector::iterator restr = restrs.begin();
    FieldRestrictionVector::iterator last_restriction = restrs.end();

    restr = std::lower_bound(restr,last_restriction,tmp);

    const bool new_restriction = ( ( restr == last_restriction ) || !(*restr == tmp) );

    if ( new_restriction ) {
      // New field restriction, verify we are not committed:
      ThrowRequireMsg(!m_meta_data->is_commit(), "mesh MetaData has been committed.");
      unsigned num_subsets = 0;
      for(FieldRestrictionVector::iterator i=restrs.begin(), iend=restrs.end(); i!=iend; ++i) {
        if (i->entity_rank() != arg_entity_rank) continue;

        const Part& partI = *m_meta_data->get_parts()[i->part_ordinal()];
        bool found_subset = contain(arg_part.subsets(), partI);
        if (found_subset) {
          ThrowErrorMsgIf( i->not_equal_stride(tmp),
            arg_method << " FAILED for " << *this << " " <<
            print_restriction( *i, arg_entity_rank, arg_part, m_field_rank ) <<
            " WITH INCOMPATIBLE REDECLARATION " <<
            print_restriction( tmp, arg_entity_rank, arg_part, m_field_rank ));
          *i = tmp;
          ++num_subsets;
        }

        bool found_superset = contain(arg_part.supersets(), partI);
        if (found_superset) {
          ThrowErrorMsgIf( i->not_equal_stride(tmp),
            arg_method << " FAILED for " << *this << " " <<
            print_restriction( *i, arg_entity_rank, arg_part, m_field_rank ) <<
            " WITH INCOMPATIBLE REDECLARATION " <<
            print_restriction( tmp, arg_entity_rank, arg_part, m_field_rank ));
          //if there's already a restriction for a superset of this part, then 
          //there's nothing to do and we're out of here..
          return;
        }
      }
      if (num_subsets == 0) {
        restrs.insert( restr , tmp );
      }
      else {
        //if subsets were found, we replaced them with the new restriction. so now we need
        //to sort and unique the vector, and trim it to remove any duplicates:
        std::sort(restrs.begin(), restrs.end());
        FieldRestrictionVector::iterator it = std::unique(restrs.begin(), restrs.end());
        restrs.resize(it - restrs.begin());
      }
    }
    else {
      ThrowErrorMsgIf( restr->not_equal_stride(tmp),
          arg_method << " FAILED for " << *this << " " <<
          print_restriction( *restr, arg_entity_rank, arg_part, m_field_rank ) <<
          " WITH INCOMPATIBLE REDECLARATION " <<
          print_restriction( tmp, arg_entity_rank, arg_part, m_field_rank ));
    }
  }
}

void FieldBaseImpl::insert_restriction(
  const char     * arg_method ,
  EntityRank       arg_entity_rank ,
  const Selector & arg_selector ,
  const unsigned * arg_stride,
  const void*      arg_init_value )
{
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::insert_restriction", LOG_FIELD, m_ordinal);

  FieldRestriction tmp( arg_entity_rank , arg_selector );

  {
    unsigned i = 0 ;
    if ( m_field_rank ) {
      for ( i = 0 ; i < m_field_rank ; ++i ) { tmp.stride(i) = arg_stride[i] ; }
    }
    else { // Scalar field is 0 == m_field_rank
      i = 1 ;
      tmp.stride(0) = 1 ;
    }
    // Remaining dimensions are 1, no change to stride
    for ( ; i < MaximumFieldDimension ; ++i ) {
      tmp.stride(i) = tmp.stride(i-1) ;
    }

    for ( i = 1 ; i < m_field_rank ; ++i ) {
      const bool bad_stride = 0 == tmp.stride(i) ||
                              0 != tmp.stride(i) % tmp.stride(i-1);
      ThrowErrorMsgIf( bad_stride,
          arg_method << " FAILED for " << *this <<
          " WITH BAD STRIDE!");
    }
  }

  if (arg_init_value != NULL) {
    //insert_restriction can be called multiple times for the same field, giving
    //the field different lengths on different mesh-parts.
    //We will only store one initial-value array, we need to store the one with
    //maximum length for this field so that it can be used to initialize data
    //for all field-restrictions. For the parts on which the field is shorter,
    //a subset of the initial-value array will be used.
    //
    //We want to end up storing the longest arg_init_value array for this field.
    //
    //Thus, we call set_initial_value only if the current length is longer
    //than what's already been stored.

    //length in bytes is num-scalars X sizeof-scalar:

    size_t num_scalars = 1;
    //if rank > 0, then field is not a scalar field, so num-scalars is
    //obtained from the stride array:
    if (m_field_rank > 0) num_scalars = tmp.stride(m_field_rank-1);

    size_t sizeof_scalar = m_data_traits.size_of;
    size_t nbytes = sizeof_scalar * num_scalars;

    size_t old_nbytes = 0;
    if (get_initial_value() != NULL) {
      old_nbytes = get_initial_value_num_bytes();
    }

    if (nbytes > old_nbytes) {
      set_initial_value(arg_init_value, num_scalars, nbytes);
    }
  }

  {
    FieldRestrictionVector & srvec = selector_restrictions();

    bool restriction_already_exists = false;
    for(FieldRestrictionVector::const_iterator it=srvec.begin(), it_end=srvec.end();
        it!=it_end; ++it) {
      if (tmp == *it) {
        restriction_already_exists = true;
        if (tmp.not_equal_stride(*it)) {
          ThrowErrorMsg("Incompatible selector field-restrictions!");
        }
      }
    }

    if ( !restriction_already_exists ) {
      // New field restriction, verify we are not committed:
      ThrowRequireMsg(!m_meta_data->is_commit(), "mesh MetaData has been committed.");
      srvec.push_back( tmp );
    }
  }
}

void FieldBaseImpl::verify_and_clean_restrictions(
  const char       * arg_method ,
  const Part& superset,
  const Part& subset,
  const PartVector & arg_all_parts )
{
  TraceIfWatching("stk::mesh::impl::FieldBaseImpl::verify_and_clean_restrictions", LOG_FIELD, m_ordinal);

  FieldRestrictionVector & restrs = restrictions();

  //Check whether both 'superset' and 'subset' are in this field's restrictions.
  //If they are, make sure they are compatible and remove the subset restriction.
  FieldRestrictionVector::iterator superset_restriction = restrs.end();
  FieldRestrictionVector::iterator subset_restriction = restrs.end();
  for (FieldRestrictionVector::iterator i = restrs.begin() ; i != restrs.end() ; ++i ) {
    if (i->part_ordinal() == superset.mesh_meta_data_ordinal()) {
      superset_restriction = i;
      if (subset_restriction != restrs.end() && subset_restriction->entity_rank() == superset_restriction->entity_rank()) break;
    }
    if (i->part_ordinal() == subset.mesh_meta_data_ordinal()) {
      subset_restriction = i;
      if (superset_restriction != restrs.end() && subset_restriction->entity_rank() == superset_restriction->entity_rank()) break;
    }
  }

  if (superset_restriction != restrs.end() && subset_restriction != restrs.end() &&
      superset_restriction->entity_rank() == subset_restriction->entity_rank()) {
    ThrowErrorMsgIf( superset_restriction->not_equal_stride(*subset_restriction),
      "Incompatible field restrictions for parts "<<superset.name()<<" and "<<subset.name());

    restrs.erase(subset_restriction);
  }
}

const void* FieldBaseImpl::get_initial_value() const
{
  return m_field_states[0]->m_impl.m_initial_value;
}

void* FieldBaseImpl::get_initial_value() {
  return m_field_states[0]->m_impl.m_initial_value;
}

unsigned FieldBaseImpl::get_initial_value_num_bytes() const {
  return m_field_states[0]->m_impl.m_initial_value_num_bytes;
}

void FieldBaseImpl::set_initial_value(const void* new_initial_value, unsigned num_scalars, unsigned num_bytes) {
  void*& init_val = m_field_states[0]->m_impl.m_initial_value;

  delete [] reinterpret_cast<char*>(init_val);
  init_val = new char[num_bytes];

  m_field_states[0]->m_impl.m_initial_value_num_bytes = num_bytes;

  m_data_traits.copy(init_val, new_initial_value, num_scalars);
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
  //(putting static here helps performance significantly but is NOT THREAD SAFE !!!)
  static FieldRestriction restr;
  restr.set_entity_rank( entity_rank );
  restr.set_part_ordinal( part.mesh_meta_data_ordinal() );

  while ( ie == ( i = find( rMap , restr ) ) && ipe != ip ) {
    // Not found try another superset part:
    restr.set_entity_rank( entity_rank );
    restr.set_part_ordinal( (*ip)->mesh_meta_data_ordinal() );
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
      const unsigned len = m_field_rank ? i->stride( m_field_rank - 1 ) : 1 ;
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

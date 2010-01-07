/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string.h>

#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

const char * field_state_name( FieldState s )
{
  static const char * name_list[] = {
    "StateNew" ,
    "StateOld" ,
    "StateNM1" ,
    "StateNM2" ,
    "StateNM3" ,
    "StateNM4" ,
    "ERROR" };

  unsigned i = s ;
  if ( StateNM4 < i ) { i = MaximumFieldStates ; }
  return name_list[i] ;
}

//----------------------------------------------------------------------

namespace {

struct RestrictionLess {
  bool operator()( const FieldBase::Restriction & lhs ,
                   const FieldBase::Restriction & rhs ) const
    { return lhs.key < rhs.key ; }

  bool operator()( const FieldBase::Restriction & lhs ,
                   const EntityKey & rhs ) const
    { return lhs.key < rhs ; }
};

std::vector<FieldBase::Restriction>::const_iterator
find( const std::vector<FieldBase::Restriction> & v ,
      const EntityKey & key )
{
  std::vector<FieldBase::Restriction>::const_iterator
    i = std::lower_bound( v.begin() , v.end() , key , RestrictionLess() );

  if ( i != v.end() && i->key != key ) { i = v.end(); }

  return i ;
}

}

//----------------------------------------------------------------------

FieldBase::~Field()
{ }

FieldBase::Field(
  MetaData * arg_mesh_meta_data ,
  unsigned   arg_ordinal ,
  const std::string & arg_name ,
  const DataTraits & arg_traits ,
  unsigned   arg_number_of_states ,
  FieldState arg_this_state )
: m_name( arg_name ),
  m_attribute(),
  m_data_traits( arg_traits ),
  m_mesh_meta_data( arg_mesh_meta_data ),
  m_mesh_meta_data_ordinal( arg_ordinal ),
  m_num_states( arg_number_of_states ),
  m_this_state( arg_this_state ),
  m_rank( 0 ),
  m_dim_map()
{
  FieldBase * const pzero = NULL ;
  const shards::ArrayDimTag * const dzero = NULL ;
  Copy<MaximumFieldStates>(    m_field_states , pzero );
  Copy<MaximumFieldDimension>( m_dim_tags ,     dzero );
}

const std::vector<FieldBase::Restriction> & FieldBase::restrictions() const
{ return m_field_states[0]->m_dim_map ; }

std::vector<FieldBase::Restriction> & FieldBase::restrictions()
{ return m_field_states[0]->m_dim_map ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void
print_field_type( std::ostream                      & arg_msg ,
                  const DataTraits                  & arg_traits ,
                  unsigned                            arg_rank ,
                  const shards::ArrayDimTag * const * arg_tags )
{
  arg_msg << "Field<" ;
  arg_msg << arg_traits.name ;
  for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
    arg_msg << "," << arg_tags[i]->name();
  }
  arg_msg << ">" ;
}

namespace {

void throw_name_suffix( const char * method ,
                        const std::string & name ,
                        const char * suffix )
{
  std::ostringstream msg ;
  msg << method ;
  msg << " FAILED: name = \"" ;
  msg << name ;
  msg << "\" CANNOT HAVE THE RESERVED STATE SUFFIX \"" ;
  msg << suffix ;
  msg << "\"" ;
  throw std::runtime_error( msg.str() );
}

// Check for compatibility:
// 1) Scalar type must match
// 2) Number of states must match
// 3) Dimension must be different by at most one rank,
//    where the tags match for the smaller rank.
void verify_field_type( const char                        * arg_method ,
                        const FieldBase                   & arg_field ,
                        const DataTraits                  & arg_traits ,
                        unsigned                            arg_rank ,
                        const shards::ArrayDimTag * const * arg_dim_tags ,
                        unsigned                    arg_num_states )
{

  const bool ok_traits = arg_traits.is_void
                      || & arg_traits == & arg_field.data_traits();

  const bool ok_number_states =
    ! arg_num_states || arg_num_states == arg_field.number_of_states();

  bool ok_dimension = ! arg_rank || arg_rank     == arg_field.rank() ||
                                    arg_rank + 1 == arg_field.rank() ||
                                    arg_rank - 1 == arg_field.rank() ;

  const unsigned check_rank = arg_rank < arg_field.rank() ?
                              arg_rank : arg_field.rank() ;

  for ( unsigned i = 0 ; i < check_rank && ok_dimension ; ++i ) {
    ok_dimension = arg_dim_tags[i] == arg_field.dimension_tags()[i] ;
  }

  if ( ! ok_traits || ! ok_number_states || ! ok_dimension ) {

    std::ostringstream msg ;             

    msg << arg_method << " FAILED: Existing field = " ;

    print_field_type( msg , arg_field.data_traits() ,
                            arg_field.rank() , 
                            arg_field.dimension_tags() ); 
 
    msg << "[ name = \"" << arg_field.name();
    msg << "\" , #states = " << arg_field.number_of_states() << " ]" ;
    msg << " Expected field info = " ;
                                           
    print_field_type( msg , arg_traits , arg_rank , arg_dim_tags );
    msg << "[ #states = " << arg_num_states << " ]" ;
   
    throw std::runtime_error( msg.str() );   
  }
}

}

FieldBase * get_field(
  const char                        * arg_method ,
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ,
  const std::vector<FieldBase*>     & arg_meta_data_fields )
{
  FieldBase * f = NULL ;

  for ( std::vector<FieldBase*>::const_iterator
        j =  arg_meta_data_fields.begin() ; 
        j != arg_meta_data_fields.end() && NULL == f ; ++j ) {
    if ( equal_case( (*j)->name() , arg_name ) ) {

      f = *j ;

      verify_field_type( arg_method , *f , arg_traits ,
                         arg_rank , arg_dim_tags , arg_num_states );
    }
  }
  return f ;
}


//----------------------------------------------------------------------

FieldBase *
FieldBase::declare_field(
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ,
  MetaData                          * arg_meta_data ,
  std::vector<FieldBase*>           & arg_meta_data_fields )
{
  static const char method[] = "stk::mesh::FieldBase::declare_field" ;

  static const char reserved_state_suffix[6][8] = {
    "_OLD" , "_N" , "_NM1" , "_NM2" , "_NM3" , "_NM4" };

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len_name   = arg_name.size();
    const int len_suffix = strlen( reserved_state_suffix[i] );
    const int offset     = len_name - len_suffix ;
    if ( 0 <= offset ) {
      const char * const name_suffix = arg_name.c_str() + offset ;
      if ( equal_case( name_suffix , reserved_state_suffix[i] ) ) {
        throw_name_suffix( method , arg_name , reserved_state_suffix[i] );
      }
    }
  }

  // Check that the field of this name has not already been declared

  FieldBase * f[ MaximumFieldStates ] ;

  f[0] = get_field( method , arg_name , arg_traits ,
                    arg_rank , arg_dim_tags , arg_num_states ,
                    arg_meta_data_fields );

  if ( NULL != f[0] ) {
    for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
      f[i] = f[0]->m_field_states[i] ; 
    }
  }
  else {
    // Field does not exist then create it

    std::string field_names[ MaximumFieldStates ];

    field_names[0] = arg_name ;

    if ( 2 == arg_num_states ) {
      field_names[1] = arg_name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
        field_names[i] = arg_name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {

      f[i] = new FieldBase( arg_meta_data , arg_meta_data_fields.size() ,
                            field_names[i] , arg_traits ,
                            arg_num_states , static_cast<FieldState>(i) );

      arg_meta_data_fields.push_back( f[i] );
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      for ( unsigned j = 0 ; j < arg_num_states ; ++j ) {
        f[i]->m_field_states[j] = f[j] ;
      }
    }
  }

  // Update the rank and dimension tags

  const unsigned old_rank = f[0]->m_rank ;

  if ( old_rank < arg_rank ) {
    for ( unsigned j = 0 ; j < arg_num_states ; ++j ) {
      f[j]->m_rank = arg_rank ;
      for ( unsigned i = old_rank ; i < arg_rank ; ++i ) {
        f[j]->m_dim_tags[i] = arg_dim_tags[i] ;
      }
    }
  }

  return f[0] ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void print_restriction( std::ostream & os ,
                        unsigned type ,
                        const Part & part ,
                        unsigned rank ,
                        const FieldBase::Restriction::size_type * stride )
{
  os << "{ entity_type(" << type << ") part(" << part.name() << ") : " ;
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

}

// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void FieldBase::insert_restriction(
  const char     * arg_method ,
  EntityType         arg_entity_type ,
  const Part     & arg_part ,
  const unsigned * arg_stride )
{
  FieldBase::Restriction tmp ;

  tmp.key = EntityKey( arg_entity_type , arg_part.mesh_meta_data_ordinal() );

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
      if ( 0 == tmp.stride[i] || 0 != tmp.stride[i] % tmp.stride[i-1] ) {
        std::ostringstream msg ;
        msg << arg_method << " FAILED for " << *this ;
        msg << " WITH BAD STRIDE " ;
        print_restriction( msg, arg_entity_type, arg_part, m_rank, tmp.stride);
        throw std::runtime_error( msg.str() );
      }
    }
  }

  {
    FieldBase::RestrictionVector & rMap = restrictions();

    FieldBase::RestrictionVector::iterator i = rMap.begin(), j = rMap.end();

    i = std::lower_bound(i,j,tmp,RestrictionLess());

    if ( i == j || i->key != tmp.key ) {
      rMap.insert( i , tmp );
    }
    else if ( Compare<MaximumFieldDimension>::not_equal(i->stride,tmp.stride) ){
      std::ostringstream msg ;
      msg << arg_method << " FAILED for " << *this << " " ;
      print_restriction( msg, arg_entity_type, arg_part, m_rank, i->stride );
      msg << " WITH INCOMPATIBLE REDECLARATION " ;
      print_restriction( msg, arg_entity_type, arg_part, m_rank, tmp.stride );
      throw std::runtime_error( msg.str() );
    }
  }
}

void FieldBase::verify_and_clean_restrictions(
  const char       * arg_method ,
  const PartVector & arg_all_parts )
{
  RestrictionVector & rMap = restrictions();
  RestrictionVector::iterator i , j ;

  for ( i = rMap.begin() ; i != rMap.end() ; ++i ) {
    if ( i->key != EntityKey() ) {
      const unsigned typeI = entity_type( i->key );
      const Part   & partI = * arg_all_parts[ entity_id( i->key ) ];
      bool  found_superset = false ;

      for ( j = i + 1 ; j != rMap.end() && ! found_superset ; ++j ) {
        if ( j->key != EntityKey() ) {
          const unsigned typeJ = entity_type( j->key );
          const Part   & partJ = * arg_all_parts[ entity_id( j->key ) ];

          if ( typeI == typeJ ) {
            const bool found_subset = contain( partI.subsets() , partJ );
            found_superset = ! found_subset &&
                             contain( partI.supersets() , partJ );

            if ( found_subset || found_superset ) {
              if ( Compare< MaximumFieldDimension >::not_equal( i->stride ,
                                                                j->stride ) ) {
                std::ostringstream msg ;
                msg << arg_method << "[" ;
                msg << *this ;
                msg << "] FAILED: " ;
                print_restriction( msg, typeI, partI, m_rank, i->stride );
                if ( found_subset ) { msg << " INCOMPATIBLE SUBSET " ; }
                else                { msg << " INCOMPATIBLE SUPERSET " ; }
                print_restriction( msg, typeJ, partJ, m_rank, j->stride );
                throw std::runtime_error( msg.str() );
              }
            }

            if ( found_subset ) { j->key = EntityKey(); }
          }
        }
        if ( found_superset ) { i->key = EntityKey(); }
      }
    }
  }

  // Clean out redundant entries:

  for ( j = i = rMap.begin() ; j != rMap.end() ; ++j ) {
    if ( j->key != EntityKey() ) {
      if ( i->key == EntityKey() ) {
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

const FieldBase::Restriction &
FieldBase::restriction( unsigned etype , const Part & part ) const
{
  static const FieldBase::Restriction empty ;

  const std::vector<FieldBase::Restriction> & rMap = restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator ie = rMap.end() ;
        std::vector<FieldBase::Restriction>::const_iterator i ;

  const PartVector::const_iterator ipe = part.supersets().end();
        PartVector::const_iterator ip  = part.supersets().begin() ;

  // Start with this part:
  EntityKey key = EntityKey( etype , part.mesh_meta_data_ordinal() );

  while ( ie == ( i = find( rMap , key ) ) && ipe != ip ) {
    // Not found try another superset part:
    key = EntityKey( etype , (*ip)->mesh_meta_data_ordinal() );
    ++ip ;
  }

  return ie == i ? empty : *i ;
}

unsigned FieldBase::max_size( unsigned entity_type ) const
{
  unsigned max = 0 ;

  const std::vector<FieldBase::Restriction> & rMap = restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator ie= rMap.end();
        std::vector<FieldBase::Restriction>::const_iterator i = rMap.begin();

  for ( ; i != ie ; ++i ) {
    if ( i->type() == entity_type ) {
      const unsigned len = m_rank ? i->stride[ m_rank - 1 ] : 1 ;
      if ( max < len ) { max = len ; }
    }
  }

  return max ;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  print_field_type( s , field.data_traits() ,
                        field.rank() , field.dimension_tags() );
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
  const PartVector & all_parts = field.mesh_meta_data().get_parts();
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();
  s << field ;
  s << " {" ;
  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << std::endl << b << "  " ;
    print_restriction( s, entity_type( i->key ),
                       * all_parts[ entity_id( i->key ) ], 
                       field.rank(), i->stride);
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


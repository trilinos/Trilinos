/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstring>
#include <iostream>
#include <sstream>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

FieldBase::~FieldBase()
{}

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

}

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  s << "FieldBase<" ;
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
  const PartVector & all_parts = field.mesh_meta_data().get_parts();
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();
  s << field ;
  s << " {" ;
  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << std::endl << b << "  " ;
    print_restriction( s, entity_rank( i->key ),
                       * all_parts[ entity_id( i->key ) ],
                       field.rank(), i->stride);
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


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
  const PartVector & all_parts = MetaData::get(field).get_parts();
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();
  s << field ;
  s << " {" ;
  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << std::endl << b << "  " ;
    i->print( s, i->rank(), * all_parts[ i->ordinal() ], field.rank() );
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


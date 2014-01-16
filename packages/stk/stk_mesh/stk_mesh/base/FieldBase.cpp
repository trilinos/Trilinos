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

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  s << "Field<" ;
  s << field.data_traits().name ;
  for ( unsigned i = 0 ; i < field.field_array_rank() ; ++i ) {
    s << "," << field.dimension_tags()[i]->name();
  }
  s << ">" ;

  s << "[ name: \"" ;
  s << field.name() ;
  s << "\" , #states: " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  s << b << field << std::endl;
  std::string indent = b;
  indent += "  ";
  print_restrictions(s, indent.c_str(), field);
  return s ;
}

std::ostream & print_restrictions(std::ostream & s ,
                                  const char * const b ,
                                  const FieldBase & field )
{
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();

  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << b;
    i->print( s, i->selector(), field.field_array_rank() );
    s << std::endl;
  }
  return s;
}

void FieldBase::set_mesh(stk::mesh::BulkData* bulk)
{
  if (m_mesh == NULL) {
    m_mesh = bulk;
  }
  else {
    ThrowRequireMsg(bulk == m_mesh, "Internal Error: Trying to use field " << name() << " on more than one bulk data");
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


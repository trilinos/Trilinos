/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/FieldState.hpp>

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

} //namespace mesh
} //namespace stk

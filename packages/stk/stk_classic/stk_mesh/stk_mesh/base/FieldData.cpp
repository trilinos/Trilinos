/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stk_mesh/base/FieldData.hpp>

namespace stk_classic {
namespace mesh {

//----------------------------------------------------------------------

const EntityDimension & EntityDimension::tag()
{ static const EntityDimension self ; return self ; }

const char * EntityDimension::name() const
{ static const char n[] = "EntityDimension" ; return n ; }

//----------------------------------------------------------------------

}//namespace mesh
}//namespace stk_classic


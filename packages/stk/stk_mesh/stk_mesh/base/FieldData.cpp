/**
 * @author H. Carter Edwards
 */

#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

const EntityDimension & EntityDimension::tag()
{ static const EntityDimension self ; return self ; }

const char * EntityDimension::name() const
{ static const char n[] = "EntityDimension" ; return n ; }

//----------------------------------------------------------------------

}//namespace mesh
}//namespace stk


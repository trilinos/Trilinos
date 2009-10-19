/**
 * @author H. Carter Edwards
 */

#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk {
namespace mesh {

const char * ElementNode::name() const
{ static const char n[] = "ElementNode" ; return n ; }

const ElementNode & ElementNode::tag()
{ static const ElementNode self ; return self ; }

const char * QuadratureTag::name() const
{ static const char n[] = "Quadrature" ; return n ; }

const QuadratureTag & QuadratureTag::tag()
{ static const QuadratureTag self ; return self ; }

const char * BasisTag::name() const
{ static const char n[] = "Basis" ; return n ; }

const BasisTag & BasisTag::tag()
{ static const BasisTag self ; return self ; }

} // namespace mesh
} // namespace stk


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

#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk_classic {
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
} // namespace stk_classic


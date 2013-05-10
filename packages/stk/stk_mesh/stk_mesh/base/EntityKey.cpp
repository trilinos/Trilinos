/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <ostream>
#include <stk_mesh/base/EntityKey.hpp>

namespace stk { namespace mesh {


std::ostream & operator << ( std::ostream & out, EntityKey key)
{
  out << "(" << static_cast<stk::topology::rank_t>(key.rank())
      << "," << key.id() << ")";

 return out;
}

}} // namespace stk::mesh


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_diag_EntityKey_hpp
#define stk_mesh_diag_EntityKey_hpp

#include <stk_util/diag/Writer.hpp>
#include <stk_mesh/base/EntityKey.hpp>


namespace stk_classic {
namespace mesh {

inline
stk_classic::diag::Writer &operator<<(stk_classic::diag::Writer &dout, const EntityKey &entity_key)  {
  return dout << entity_rank(entity_key) << ":" << entity_id(entity_key);
}

} // namespace stk_classic
} // namespace mesh

#endif // stk_mesh_diag_EntityKey_hpp

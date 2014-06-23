/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_BASE_BULK_MODIFICATION_HPP
#define STK_MESH_BASE_BULK_MODIFICATION_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace stk_classic {
namespace mesh {

  /** \brief Determine closure of the entities vector
   *
   *   \param bulk  BulkData must be in a parallel consistent state.
   *
   *   \param entities Each entity must be in the locally_used part.
   *
   *   \param entities_closure Parallel consistent closure of the input
   *     vector. This vector will be sorted and unique.  May include
   *     ghosted entities.
   */
void find_closure( const BulkData & bulk,
    const std::vector< Entity *> & entities,
    std::vector< Entity *> & entities_closure);

} // namespace mesh
} // namespace stk_classic


#endif // STK_MESH_BASE_BULK_MODIFICATION_HPP

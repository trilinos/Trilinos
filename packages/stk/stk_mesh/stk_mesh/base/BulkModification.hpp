#ifndef STK_MESH_BASE_BULK_MODIFICATION_HPP
#define STK_MESH_BASE_BULK_MODIFICATION_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

  /** \brief Determine closure of the entities vector
   *
   *   \param bulk  BulkData must be in a parallel consistent state.
   *
   *   \param entities Each entity must be in the locally_used part.
   *
   *   \param entities_closure Parallel consisent closure of the input
   *     vector. This vector will be sorted and unique.  May include
   *     ghosted entities.
   */
void find_closure( const BulkData & bulk,
    const std::vector< Entity *> & entities,
    std::vector< Entity *> & entities_closure);

} // namespace mesh
} // namespace stk


#endif // STK_MESH_BASE_BULK_MODIFICATION_HPP

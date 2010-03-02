/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_GetEntities_hpp
#define stk_mesh_GetEntities_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------

/** \brief  Local count selected entities of each type.
 *
 * \param selector
 * \param mesh
 * \param count  
 */
void count_entities(
  const Selector & selector ,
  const BulkData & mesh ,
  std::vector<EntityType> & count );

/** \brief Get all entities of the specified type, sorted by ID.  */
void get_entities( const BulkData & , EntityType type,
                   std::vector< Entity*> & );

/** \brief  Count entities in selected parts (selected by the
 *          given selector instance), and sorted by ID.
 */
unsigned count_selected_entities(
    const Selector & selector ,
    const std::vector< Bucket * > & input_buckets );

/** \brief  Get entities in selected parts (selected by the
 *          given selector instance), and sorted by ID.
 */
void get_selected_entities(
    const Selector & selector ,
    const std::vector< Bucket * > & input_buckets ,
    std::vector< Entity * > & entities );

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif //  stk_mesh_GetEntities_hpp

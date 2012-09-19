/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_CreateAdjacentEntities_hpp
#define stk_mesh_CreateAdjacentEntities_hpp

#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

/** Create all the internal entities and their relations,
 * i.e. given a mesh with element and node relations
 *  * create all the side and edges
 *  * and setup the relations
 *
 * Expensive function that should only be called once
 *
 * This is a parallel collective function (it should be called on all
 * processors at the same time
 *
 * \param mesh A consisent mesh with Element and node relations
 * \param add_parts Newly created entities will be added to the add_parts and
 *                  the root part for the given subcell topology.  Note, no part in the
 *                  'add_parts' vector can have a primary_entity_rank
 */
void create_adjacent_entities(
    BulkData & mesh,
    PartVector & add_parts
    );

}
}
#endif

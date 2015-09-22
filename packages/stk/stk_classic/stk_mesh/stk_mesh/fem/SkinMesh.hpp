/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_SkinMesh_hpp
#define stk_mesh_SkinMesh_hpp

#include <vector>

namespace stk_classic {
namespace mesh {

class BulkData;
class Part;
class Entity;

typedef std::vector<Entity *> EntityVector;

/**
 * Skin the entire mesh.
 */
void skin_mesh( BulkData & mesh,
                EntityRank entity_rank,
                Part * skin_part = NULL );

/**
 * Given a vector of modified/created elements, update the skin.
 */
void reskin_mesh( BulkData & mesh,
                  EntityRank entity_rank,
                  EntityVector & owned_modified_elements,
                  Part * skin_part = NULL );



}
}
#endif

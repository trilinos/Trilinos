/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_fem_EntityRanksEnums_hpp
#define stk_mesh_fem_EntityRanksEnums_hpp

#include <string>
#include <vector>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */
/** \brief  Enumeration of types of entities.
 *
 *  The types of entities is intended to be modifiable / extensible by
 *  - maintaining the valid types in a contiguous span [0..EntityRankEnd),
 *  - having the first four values correspond to the topological
 *    entity types.
 */
// DEPRECATED: 09/15/10 FEM refactor
enum EntityRankEnum {
  Node                = 0 ,
  Edge                = 1 ,
  Face                = 2 ,
  Element             = 3 ,
  Particle            = 4 ,
  Constraint          = 5 ,
  EntityRankEnd       = 6 ,
  EntityRankUndefined = -1
};

/** \brief  Finite element entity-type names */
// DEPRECATED: 09/15/10 FEM refactor
const std::vector<std::string> & fem_entity_rank_names();

// DEPRECATED: 09/15/10 FEM refactor
inline
EntityRankEnum fem_entity_rank( int t )
{ return 0 <= t && t < EntityRankEnd ? EntityRankEnum(t) : EntityRankUndefined ; }

/** \} */

}//namespace mesh
}//namespace stk

#endif


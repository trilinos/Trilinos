/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_fem_EntityTypesEnums_hpp
#define stk_mesh_fem_EntityTypesEnums_hpp

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
 *  - maintaining the valid types in a contiguous span [0..EntityTypeEnd),
 *  - having the first four values correspond to the topological
 *    entity types.
 */
enum EntityTypeEnum {
  Node                = 0 ,
  Edge                = 1 ,
  Face                = 2 ,
  Element             = 3 ,
  Particle            = 4 ,
  Constraint          = 5 ,
  EntityTypeEnd       = 6 ,
  EntityTypeUndefined = -1
};

/** \brief  Finite element entity-type names */
const std::vector<std::string> & fem_entity_type_names();

inline
EntityTypeEnum fem_entity_type( int t )
{ return 0 <= t && t < EntityTypeEnd ? EntityTypeEnum(t) : EntityTypeUndefined ; }

/** \} */

}//namespace mesh
}//namespace stk

#endif


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_BoundaryAnalysis_hpp
#define stk_mesh_BoundaryAnalysis_hpp

#include <vector>

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/fem/FEMInterface.hpp>

#include <stk_util/util/NamedPair.hpp>

struct CellTopologyData;

namespace stk {
namespace mesh {

class BulkData;
class Entity;

/**
 * A pair of Entity* and a local side id (defined a part of a side)
 */
NAMED_PAIR(EntitySideComponent, Entity *, entity, unsigned, side_ordinal)
/**
 * A pair of EntitySideComponents (defines a side of the boundary)
 * Most sides will have two EntitySideComponents, but degenerate cases
 * ie shells will can have sides with more than two components)
 */
NAMED_PAIR(EntitySide, EntitySideComponent, inside, EntitySideComponent, outside)

/**
 * A vector of EntitySide (defines a boundary)
 */
typedef std::vector<EntitySide> EntitySideVector;

typedef std::vector<EntitySideComponent> EntitySideComponentVector;

/** \brief Given a closure, return a boundary of items of closure_rank-1
 *
 * A boundary will contain the entities "touching" the "outside" of the
 * closure.
 */
void boundary_analysis(const BulkData & bulk_data,
                       const EntityVector & entities_closure,
                       EntityRank closure_rank,
                       EntitySideVector& boundary);


}
}
#endif

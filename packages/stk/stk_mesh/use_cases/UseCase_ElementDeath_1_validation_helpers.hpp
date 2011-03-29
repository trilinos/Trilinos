/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_use_cases_UseCase_ElementDeath_1_validation_helpers_hpp
#define stk_mesh_use_cases_UseCase_ElementDeath_1_validation_helpers_hpp

#include <stk_util/parallel/Parallel.hpp>


namespace stk {
namespace mesh {

class Entity;
class BulkData;

namespace fixtures {

class GridFixture;

}

}
}

//Generates a vector of entities to be killed in this iteration
std::vector<stk::mesh::Entity *> entities_to_be_killed( 
    const stk::mesh::BulkData & mesh, 
    int iteration, 
    stk::mesh::EntityRank entity_rank
    );

//Validates that the correct entites were killed in this iteration
bool validate_iteration(
    stk::ParallelMachine pm,
    stk::mesh::fixtures::GridFixture & fixture,
    int iteration
    );

#endif

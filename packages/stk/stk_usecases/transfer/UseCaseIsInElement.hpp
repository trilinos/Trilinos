/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <utility>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

namespace stk {
namespace usecase {

typedef stk::mesh::Field<double, stk::mesh::Cartesian>       VectorField ;

size_t is_in_element(VectorField *domain_coordinates,
                     VectorField *range_coordinates,
                     const std::vector<std::pair<stk::mesh::Entity*, stk::mesh::Entity*> > &entity_map,
                     std::vector<std::size_t> &not_in_element);
}
}

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_REFINENEARLEVELSETS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_REFINENEARLEVELSETS_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <functional>

namespace krino {

void refine_elements_that_intersect_distance_interval_from_levelset(stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const stk::mesh::Field<double> & stkLevelSetField,
    const std::function<void()> & initialize_levelsets,
    const std::array<double,2> & refinementDistanceInterval,
    const unsigned numRefinementLevels);
}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_REFINENEARLEVELSETS_HPP_ */

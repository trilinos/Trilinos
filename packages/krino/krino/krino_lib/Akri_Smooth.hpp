/*
 * Akri_Smooth.hpp
 *
 *  Created on: Nov 8, 2024
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

void improve_quality_by_ODT_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations = 50);

void improve_quality_by_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations = 50);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_ */

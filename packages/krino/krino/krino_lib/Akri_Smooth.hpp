#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_DistributedVector.hpp>

namespace krino {

using NodeSearchDirectionFilter = std::function<void(stk::mesh::Entity, stk::math::Vector3d&)>;
using NodeObjFn = std::function<double(const stk::math::Vector3d&)>;
using NodeObjFnSens = std::function<void(const stk::math::Vector3d&, stk::math::Vector3d&)>;
using MeshNodesObjFn = std::function<double(const DistributedVector&)>;
using MeshNodesObjFnSens = std::function<void(const DistributedVector&, DistributedVector&)>;

using OptimizeNodeLocation = std::function<void(const NodeObjFn &objFn, const NodeObjFnSens &gradObjFn, stk::math::Vector3d& nodeLocation)>;
using OptimizeMeshNodeLocations = std::function<void(const MeshNodesObjFn &objFn, const MeshNodesObjFnSens &gradObjFn, DistributedVector& nodeLocations)>;

void improve_quality_by_ODT_smoothing_on_interior(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const size_t maxNumSmoothIterations = 50);

void improve_quality_by_optimized_mean_ratio_smoothing_on_interior(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeNodeLocation & optimize_node_location,
    const size_t maxNumSmoothIterations = 50);

void improve_quality_by_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeNodeLocation & optimize_node_location,
    const NodeSearchDirectionFilter & node_search_direction_filter,
    const size_t maxNumSmoothIterations = 50);

void improve_quality_by_simultaneous_optimized_mean_ratio_smoothing(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizeMeshNodeLocations & optimize_node_locations,
    const NodeSearchDirectionFilter & node_search_direction_filter);
}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SMOOTH_HPP_ */

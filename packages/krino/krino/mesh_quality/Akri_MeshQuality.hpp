#ifndef KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_
#include <stk_math/StkVector.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }

namespace krino
{
class QualityMetric;
class CoordinatesFieldRef;
class FieldRef;

double compute_mesh_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const QualityMetric &qualityMetric);

template<typename ELEMCONTAINER>
double compute_minimum_element_quality(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const ELEMCONTAINER & elements,
    const QualityMetric &qualityMetric);

double compute_quality_if_node_is_moved_terminating_early_if_below_threshold(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    stk::mesh::Entity node,
    const stk::math::Vector3d & newLocation,
    const QualityMetric &qualityMetric,
    const double qualityThreshold);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_ */

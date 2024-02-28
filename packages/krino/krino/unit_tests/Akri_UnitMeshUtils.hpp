#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

class FieldRef;

stk::mesh::Entity find_local_node_closest_to_location(const stk::mesh::BulkData& mesh, const stk::mesh::Selector& nodeSelector, const FieldRef coordsField, const stk::math::Vector3d& location);

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNITMESHUTILS_HPP_ */

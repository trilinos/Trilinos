#ifndef KRINO_KRINO_KRINO_LIB_AKRI_ORIENTEDSIDENODES_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_ORIENTEDSIDENODES_HPP_
#include <stk_mesh/base/Types.hpp>
#include <array>

namespace krino {

std::array<stk::mesh::Entity,3> get_oriented_triangle_side_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side);

std::array<stk::mesh::Entity,2> get_oriented_line_side_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & negativeSideElementSelector, const stk::mesh::Entity side);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_ORIENTEDSIDENODES_HPP_ */

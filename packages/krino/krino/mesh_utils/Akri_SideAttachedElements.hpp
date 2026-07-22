#ifndef KRINO_KRINO_MESH_UTILS_AKRI_SIDEATTACHEDELEMENTS_HPP_
#define KRINO_KRINO_MESH_UTILS_AKRI_SIDEATTACHEDELEMENTS_HPP_
#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace krino {

std::vector<bool> are_elements_side_attached_to_selected_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector);

std::vector<stk::mesh::Entity> get_selected_owned_side_unattached_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector);

std::vector<stk::mesh::Entity> find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector);

}



#endif /* KRINO_KRINO_MESH_UTILS_AKRI_SIDEATTACHEDELEMENTS_HPP_ */

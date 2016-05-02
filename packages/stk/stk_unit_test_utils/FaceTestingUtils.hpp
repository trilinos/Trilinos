#ifndef FACETESTINGUTILS_HPP_
#define FACETESTINGUTILS_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <string>
#include <vector>

unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_count_sides(std::string filename);

unsigned read_file_count_sides(std::string filename);

bool fully_connected_elements_to_faces(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_fully_connected_stk(std::string filename);

unsigned read_file_fully_connected_stk(std::string filename);

unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename);

unsigned read_file_shared_faces_different_elements_stk(std::string filename);

unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename);

unsigned read_file_shared_faces_same_elements_stk(std::string filename);

bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::set<unsigned>& counts);

bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

namespace stk
{
namespace unit_test_util
{
stk::mesh::Entity declare_element_to_sub_topology_with_nodes(stk::mesh::BulkData &mesh, stk::mesh::Entity elem, const stk::mesh::EntityVector &sub_topology_nodes,
                        stk::mesh::EntityId global_sub_topology_id, stk::mesh::EntityRank to_rank, stk::mesh::Part &part);
}
}

#endif // FACETESTINGUTILS_HPP_

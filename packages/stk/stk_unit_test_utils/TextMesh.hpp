#ifndef TextMesh_hpp
#define TextMesh_hpp

#include <string>
#include <vector>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace unit_test_util
{
typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinatesField;
void fill_mesh_using_text_mesh(const std::string &meshDesc, stk::mesh::BulkData &bulkData);
void fill_mesh_using_text_mesh_with_coordinates(const std::string &meshDesc, const std::vector<double> &coordinates, stk::mesh::BulkData &bulkData);

} // namespace unit_test_util
} // namespace stk


#endif

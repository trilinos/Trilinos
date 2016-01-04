#ifndef unittestMeshUtils_hpp
#define unittestMeshUtils_hpp

#include <stk_mesh/base/Types.hpp>      // for EntityVector, PartVector
#include <string>                       // for string
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk
{

namespace unit_test_util
{

void put_mesh_into_part(stk::mesh::BulkData& bulkData, stk::mesh::Part& part);

std::string get_name_of_generated_mesh(int xdim, int ydim, int zdim, const std::string &options);

void move_killed_elements_out_of_parts(stk::mesh::BulkData& bulkData,
                                  const stk::mesh::EntityVector& killedElements,
                                  const stk::mesh::PartVector& removeParts);

} // namespace unit_test_util
} // namespace stk


#endif // unittestMeshUtils_hpp

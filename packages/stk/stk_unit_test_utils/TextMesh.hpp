#ifndef TextMesh_hpp
#define TextMesh_hpp

#include <string>
namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace unit_test_util
{

void fill_mesh_using_text_mesh(const std::string &meshDesc, stk::mesh::BulkData &bulkData);

} // namespace unit_test_util
} // namespace stk


#endif

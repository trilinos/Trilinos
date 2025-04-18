#ifndef STK_TOOLS_MESH_TOOLS_FIXNODESHARINGVIASEARCH_HPP_
#define STK_TOOLS_MESH_TOOLS_FIXNODESHARINGVIASEARCH_HPP_
namespace stk { namespace mesh { class BulkData; }}

namespace stk {
namespace tools {

void fix_node_sharing_via_search(stk::mesh::BulkData& bulkData, double tolerance = 1.e-6);

}}

#endif

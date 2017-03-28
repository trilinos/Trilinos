#include <vector>
#include <limits>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; }}

namespace stk {
namespace balance {
namespace internal {

std::vector<unsigned> assign_target_subdomains_roundrobin_to_procs(unsigned num_procs_M, unsigned num_procs_N);
void fill_decomp(const int num_partitions, stk::mesh::BulkData& bulk, stk::mesh::EntityProcVec &decomp);
stk::mesh::EntityProcVec get_element_decomp(const int num_partitions, stk::mesh::BulkData& bulk);

}}}

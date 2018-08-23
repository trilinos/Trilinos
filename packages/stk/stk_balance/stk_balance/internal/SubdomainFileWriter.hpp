#ifndef STK_STK_BALANCE_STK_BALANCE_INTERNAL_SUBDOMAINFILEWRITER_HPP_
#define STK_STK_BALANCE_STK_BALANCE_INTERNAL_SUBDOMAINFILEWRITER_HPP_
#include <stk_util/parallel/Parallel.hpp>
#include <string>
#include <tuple>

namespace stk { namespace mesh { class BulkData; }}

namespace stk {
namespace balance {
namespace internal {

std::tuple<int,int> get_included_and_num_target_procs(stk::mesh::BulkData &bulk, stk::ParallelMachine comm);

int get_subdomain_index(int includeMe, stk::ParallelMachine comm);

void write_subdomain_files(stk::mesh::BulkData &bulk, int numTarget, int mySubdomain, const std::string& outputMesh);

}}}
#endif /* STK_STK_BALANCE_STK_BALANCE_INTERNAL_SUBDOMAINFILEWRITER_HPP_ */

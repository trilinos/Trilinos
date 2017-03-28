#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityCommListInfo.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <iostream>

namespace stk {
namespace mesh {
namespace impl {

bool is_comm_list_globally_consistent(stk::ParallelMachine communicator, const EntityCommListInfoVector& comm_list);
bool is_comm_list_globally_consistent(stk::ParallelMachine communicator, const EntityCommListInfoVector& comm_list, std::ostream& error_msg);

} //namespace impl
} //namespace mesh
} //namespace stk


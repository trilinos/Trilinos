
#ifndef DIFFINGTOOL_HPP_
#define DIFFINGTOOL_HPP_

#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <string>                       // for string
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { class CommBuffer; }
namespace stk { class CommSparse; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { struct EnvData; }

namespace stk {
namespace diff {

//----- Functions for the App to call ----------------

bool parts_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data);
bool parts_match_except(const stk::mesh::BulkData& bulk, stk::EnvData& env_data, stk::mesh::Part* skipPart);


//------- Implementation functions -------------

void pack_string(stk::CommBuffer& buf, const std::string& name);
std::string unpack_string(stk::CommBuffer& buf);
void pack_part_names(stk::CommBuffer& buf, const stk::mesh::PartVector& parts);
stk::CommBuffer& get_comm_buffer_for_destination_proc(stk::CommSparse& comm);
void send_part_names_to_diffing_tool(const stk::mesh::BulkData& bulk, stk::ParallelMachine communicator);
void allocate_or_communicate(int iphase, stk::CommSparse& comm);
int get_global_part_differences_for_app(stk::ParallelMachine comm);
int get_global_part_differences(stk::ParallelMachine comm, int numLocalDiffs);
int parallel_sum(stk::ParallelMachine comm, int numLocal);
bool bucket_part_memberships_match(const stk::mesh::BulkData& bulk, stk::EnvData& env_data);
int get_global_bucket_part_membership_differences(stk::ParallelMachine comm, int numLocalDiffs);
int get_global_bucket_count_differences(stk::ParallelMachine comm, int numLocalDiffs);
void communicate_run_state(stk::EnvData& env_data, bool continue_runs);

}
}

#endif

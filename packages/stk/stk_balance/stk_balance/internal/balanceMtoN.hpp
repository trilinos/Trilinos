#ifndef STK_BALANCE_MTON_HPP
#define STK_BALANCE_MTON_HPP

#include <string>

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace balance {
namespace internal {

void rebalanceMtoN(stk::mesh::BulkData& bulkData, int num_target_procs, const std::string& outputFilename, int numSteps = -1, double timeStep = 0.0);

}}}

#endif

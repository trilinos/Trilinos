#ifndef STK_BALANCE_MTON_HPP
#define STK_BALANCE_MTON_HPP

#include <stk_mesh/base/Field.hpp>
#include <string>

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace balance {
namespace internal {

void rebalanceMtoN(stk::mesh::BulkData& bulkData, stk::mesh::Field<double> &targetDecompField, int num_target_procs, const std::string& outputFilename, int numSteps = -1, double timeStep = 0.0);

}}}

#endif

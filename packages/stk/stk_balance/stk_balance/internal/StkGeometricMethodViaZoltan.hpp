#ifndef STK_GEOMETRIC_METHOD_VIA_ZOLTAN_HPP
#define STK_GEOMETRIC_METHOD_VIA_ZOLTAN_HPP

#include <vector>
#include <string>
#include <mpi.h>

namespace stk { namespace balance { namespace internal { class GeometricVertices; }}}

namespace stk {
namespace balance {

class BalanceSettings;

std::vector<unsigned> get_decomposition(const stk::balance::internal::GeometricVertices& vertexInfo, const BalanceSettings& balanceSettings, int numParts, MPI_Comm comm);

}
}
#endif // STK_GEOMETRIC_METHOD_VIA_ZOLTAN_HPP

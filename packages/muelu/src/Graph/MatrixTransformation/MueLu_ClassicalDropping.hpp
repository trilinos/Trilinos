#ifndef MUELU_CLASSICALDROPPING_HPP
#define MUELU_CLASSICALDROPPING_HPP

#include "MueLu_DroppingCommon.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace MueLu::ClassicalDropping {

template <class local_matrix_type, class diag_view_type>
class DropFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using memory_space       = typename local_matrix_type::memory_space;
  using results_view       = Kokkos::View<DecisionType*, memory_space>;

  using ATS                 = Kokkos::ArithTraits<scalar_type>;
  using magnitudeType       = typename ATS::magnitudeType;
  using boundary_nodes_view = Kokkos::View<bool*, memory_space>;

  local_matrix_type A;
  diag_view_type diag;  // corresponds to overlapped diagonal
  magnitudeType eps;
  boundary_nodes_view boundaryNodes;
  results_view results;

 public:
  DropFunctor(local_matrix_type& A_, magnitudeType threshold, diag_view_type diag_, boundary_nodes_view boundaryNodes_, results_view& results_)
    : A(A_)
    , diag(diag_)
    , eps(threshold)
    , boundaryNodes(boundaryNodes_)
    , results(results_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(const local_ordinal_type rlid) {
    auto row      = A.rowConst(rlid);
    size_t offset = A.graph.row_map(rlid);
    for (local_ordinal_type k = 0; k < row.length; ++k) {
      auto clid = row.colidx(k);

      auto val    = row.value(k);
      auto aiiajj = ATS::magnitude(diag(rlid)) * ATS::magnitude(diag(clid));  // |a_ii|*|a_jj|
      auto aij2   = ATS::magnitude(val) * ATS::magnitude(val);                // |a_ij|^2

      results(offset + k) = (aij2 <= eps * eps * aiiajj) ? DROP : KEEP;
    }
    return false;
  }
};

}  // namespace MueLu::ClassicalDropping

#endif

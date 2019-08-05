#ifndef TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DECL_HPP
#define TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DECL_HPP

#include "Tpetra_RowMatrix_fwd.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <tuple>

namespace Tpetra {
namespace Details {

namespace Impl {
  template <class NT>
  struct LocalRowOffsetsResult {
  private:
    using local_graph_type =
      typename KokkosSparse::CrsMatrix<
        double, int, typename NT::execution_space, void>::
          staticcrsgraph_type;
  public:
    using offsets_type =
      typename local_graph_type::row_map_type::non_const_type;
    using offset_type = typename offsets_type::non_const_value_type;

    offsets_type ptr;
    offset_type nnz;
    size_t maxNumEnt;
  };
} // namespace Impl

template <class SC, class LO, class GO, class NT>
Impl::LocalRowOffsetsResult<NT>
localRowOffsets (const RowMatrix<SC, LO, GO, NT>& A);

template <class SC, class LO, class GO, class NT>
KokkosSparse::CrsMatrix<
  typename Kokkos::ArithTraits<SC>::val_type,
    LO,
    typename NT::execution_space,
    void>
localDeepCopyLocallyIndexedRowMatrix
  (const RowMatrix<SC, LO, GO, NT>& A,
   const char label[]);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DECL_HPP

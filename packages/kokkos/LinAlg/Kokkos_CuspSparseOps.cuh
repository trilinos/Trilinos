#include "Kokkos_CuspSparseOps.hpp"
// #include <cusp/multiply.h>
// #include <cusp/detail/device/spmv/hyb.h>

namespace Kokkos {

  template<class Scalar, class Ordinal, class Node>
  template <class O, class S>
  static void CuspSparseOps<Scalar,Ordinal,Node>::cusp_mult(const cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> &mat, 
                                                            const Scalar *X, size_t ldx, Scalar *Y, size_t ldy) 
  {
    assert(false);
  }

  template<class Scalar, class Ordinal, class Node>
  template <class O, class S> 
  static void CuspSparseOps<Scalar,Ordinal,Node>::cusp_convert(const cusp::crs_matrix<Ordinal,Scalar,cusp::device_memory> *crs
                              cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> &hyb)
  {
    assert(false);
  }


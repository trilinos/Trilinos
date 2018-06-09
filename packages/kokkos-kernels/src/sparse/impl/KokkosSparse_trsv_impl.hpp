/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSSPARSE_IMPL_TRSM_HPP_
#define KOKKOSSPARSE_IMPL_TRSM_HPP_

/// \file Kokkos_Sparse_impl_trsm.hpp
/// \brief Implementation(s) of sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <vector> // temporarily

namespace KokkosSparse {
namespace Impl {
namespace Sequential {

template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsrUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type r = 0; r < numRows; ++r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // for each row r
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsr (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type r = 0; r < numRows; ++r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }

    matrix_scalar_type A_rr = STS::zero ();
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);

    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  The
      // way to fix this would be for Tpetra::CrsMatrix to remember
      // the local column index of the diagonal entry (if there is
      // one) in each row, and pass that along to this function.
      if (r == c) {
        A_rr += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // for each row r
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsrUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  // If local_ordinal_type is unsigned and numRows is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numRows == 0) {
    return;
  }

  // Don't use r >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do r == 0 (last
  // iteration) below.
  for (local_ordinal_type r = numRows - 1; r != 0; --r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // for each row r

  // Last iteration: r = 0.
  {
    const local_ordinal_type r = 0;
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // last iteration: r = 0
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsr (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  // If local_ordinal_type is unsigned and numRows is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numRows == 0) {
    return;
  }

  // Don't use r >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do r == 0 (last
  // iteration) below.
  for (local_ordinal_type r = numRows - 1; r != 0; --r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    // We assume the diagonal entry is first in the row.
    const matrix_scalar_type A_rr = val(beg);
    for (offset_type k = beg + static_cast<offset_type> (1); k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // for each row r

  // Last iteration: r = 0.
  {
    const local_ordinal_type r = 0;
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    // We assume the diagonal entry is first in the row.
    const matrix_scalar_type A_rr = val(beg);
    for (offset_type k = beg + 1; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // last iteration: r = 0
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type r = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type r = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsc (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscUnitDiagConj (RangeMultiVectorType X,
                              const CrsMatrixType& A,
                              DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscConj (RangeMultiVectorType X,
                      const CrsMatrixType& A,
                      DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsc (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscUnitDiagConj (RangeMultiVectorType X,
                              const CrsMatrixType& A,
                              DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscConj (RangeMultiVectorType X,
                      const CrsMatrixType& A,
                      DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.extent(1);
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c
}

} // namespace Sequential
} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOSSPARSE_IMPL_TRSM_HPP

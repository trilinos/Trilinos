//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSSPARSE_IMPL_SOR_HPP
#define KOKKOSSPARSE_IMPL_SOR_HPP

/// \file KokkosSparse_impl_sor.hpp
/// \brief Sequential implementations of Gauss-Seidel and SOR.
///
/// This file exists mainly as a temporary porting aid.  Until we can
/// set up thread-parallel versions of these kernels, we have
/// sequential versions here.  They are not terribly well optimized.
/// The point is to have reasonable sequential defaults, not
/// super-fast kernels.  If you want super-fast kernels, reimplement
/// them in Kokkos and try to make them parallel.
///
/// The exception to this might be Gauss-Seidel.  Gauss-Seidel is
/// harder to parallelize without changing the semantics, unless we
/// implement it with triangular solves and do the latter in parallel.
/// There is a set-up cost for thread-parallel sparse triangular
/// solves, and speed-up isn't even close to perfect, so it might pay
/// for smaller thread counts to have an optimized sequential kernel.
/// We have <i>not</i> done this here.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <vector>  // temporarily

namespace KokkosSparse {
namespace Impl {
namespace Sequential {

/// \brief Implementation of local Gauss-Seidel.
///
/// "Local" means local to the MPI process.
///
/// \tparam LocalOrdinal The type of each column index.
/// \tparam OffsetType The type of the entries of the row offsets array.
/// \tparam MatrixScalar The type of the entries (values of the matrix.
/// \tparam DomainScalar The type of the entries of the input multivector.
/// \tparam RangeScalar The type of the entries of the output multivector.
///
/// \param numRows [in] Number of rows in the (local) matrix.
/// \param numCols [in] Number of columns in the input multivector B.
///   <i>NOT</i> the number of columns in the matrix!
/// \param ptr [in] The matrix's row offsets.
/// \param ind [in] The matrix's column indices.
/// \param val [in] The matrix's values.
/// \param B [in] The input multivector.
/// \param b_stride [in] Column stride of the input multivector B.
/// \param X [in] The output multivector.
/// \param x_stride [in] Column stride of the output multivector X.
/// \param D [in] Array of the "diagonal entries" of the matrix.
///   These may differ from the matrix's actual diagonal entries, as
///   in L1 Gauss-Seidel for example.
/// \param omega [in] Damping parameter.
/// \param direction [in] Sweep direction: "F" for forward, "B" for
///   backward.
template <class LocalOrdinal, class OffsetType, class MatrixScalar, class DomainScalar, class RangeScalar>
void gaussSeidel(const LocalOrdinal numRows, const LocalOrdinal numCols, const OffsetType* const ptr,
                 const LocalOrdinal* const ind, const MatrixScalar* const val, const DomainScalar* const B,
                 const OffsetType b_stride, RangeScalar* const X, const OffsetType x_stride,
                 const MatrixScalar* const D, const MatrixScalar omega, const char direction[]) {
  using Kokkos::ArithTraits;
  typedef LocalOrdinal LO;
  const OffsetType theNumRows = static_cast<OffsetType>(numRows);
  const OffsetType theNumCols = static_cast<OffsetType>(numCols);

  if (numRows == 0 || numCols == 0) {
    return;  // Nothing to do.
  } else if (numRows > 0 && ptr[numRows] == 0) {
    // All the off-diagonal entries of A are zero, and all the
    // diagonal entries are (implicitly) 1.  Therefore compute: X :=
    // (1 - omega) X + omega B.  There's no need to care about the
    // direction, since there are no cross-row data dependencies in
    // this case.
    const MatrixScalar oneMinusOmega = ArithTraits<MatrixScalar>::one() - omega;
    for (OffsetType j = 0; j < theNumCols; ++j) {
      RangeScalar* const x_j        = X + j * x_stride;
      const DomainScalar* const b_j = B + j * b_stride;
      for (OffsetType i = 0; i < theNumRows; ++i) {
        x_j[i] = oneMinusOmega * x_j[i] + omega * b_j[i];
      }
    }
    return;
  }

  if (numCols == 1) {
    if (direction[0] == 'F' || direction[0] == 'f') {
      for (LO i = 0; i < numRows; ++i) {
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    } else if (direction[0] == 'B' || direction[0] == 'b') {
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO i = numRows - 1; i != 0; --i) {
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
      {  // last loop iteration
        const LO i         = 0;
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    }
  } else {  // numCols > 1
    // mfh 20 Dec 2012: If Gauss-Seidel for multivectors with
    // multiple columns becomes important, we can add unrolled
    // implementations.  The implementation below is not unrolled.
    // It may also be reasonable to parallelize over right-hand
    // sides, if there are enough of them, especially if the matrix
    // fits in cache.
    std::vector<RangeScalar> temp(numCols);
    RangeScalar* const x_temp = numCols == 0 ? NULL : &temp[0];

    if (direction[0] == 'F' || direction[0] == 'f') {
      for (LO i = 0; i < numRows; ++i) {
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
    } else if (direction[0] == 'B' || direction[0] == 'b') {  // backward mode
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO i = numRows - 1; i != 0; --i) {
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
      {  // last loop iteration
        const LO i = 0;
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
    }
  }
}

/// \brief Implementation of reordered local Gauss-Seidel.
///
/// "Local" means local to the MPI process.
///
/// \tparam LocalOrdinal The type of each column index.
/// \tparam OffsetType The type of the entries of the row offsets array.
/// \tparam MatrixScalar The type of the entries (values of the matrix.
/// \tparam DomainScalar The type of the entries of the input multivector.
/// \tparam RangeScalar The type of the entries of the output multivector.
///
/// \param numRows [in] Number of rows in the (local) matrix.
/// \param numCols [in] Number of columns in the input multivector B.
///   <i>NOT</i> the number of columns in the matrix!
/// \param ptr [in] The matrix's row offsets.
/// \param ind [in] The matrix's column indices.
/// \param val [in] The matrix's values.
/// \param B [in] The input multivector.
/// \param b_stride [in] Column stride of the input multivector B.
/// \param X [in] The output multivector.
/// \param x_stride [in] Column stride of the output multivector X.
/// \param D [in] Array of the "diagonal entries" of the matrix.
///   These may differ from the matrix's actual diagonal entries, as
///   in L1 Gauss-Seidel for example.
/// \param rowInd [in] Array of row indices to process.  It has
///   numRowInds entries.  This array determines the order in which
///   the rows are accessed.  It is legal for this to contain fewer
///   entries than the number of rows in the matrix.
/// \param numRowInds [in] Number of entries in rowInd; the number of
///   rows to process.  This may be less than or equal to numRows (the
///   number of rows in the matrix).
/// \param omega [in] Damping parameter.
/// \param direction [in] Sweep direction: "F" for forward, "B" for
///   backward.
template <class LocalOrdinal, class OffsetType, class MatrixScalar, class DomainScalar, class RangeScalar>
void reorderedGaussSeidel(const LocalOrdinal numRows, const LocalOrdinal numCols, const OffsetType* const ptr,
                          const LocalOrdinal* const ind, const MatrixScalar* const val, const DomainScalar* const B,
                          const OffsetType b_stride, RangeScalar* const X, const OffsetType x_stride,
                          const MatrixScalar* const D, const LocalOrdinal* const rowInd,
                          const LocalOrdinal numRowInds,  // length of rowInd
                          const MatrixScalar omega, const char direction[]) {
  using Kokkos::ArithTraits;
  typedef LocalOrdinal LO;
  const OffsetType theNumRows = static_cast<OffsetType>(numRows);
  const OffsetType theNumCols = static_cast<OffsetType>(numCols);

  if (numRows == 0 || numCols == 0) {
    return;  // Nothing to do.
  } else if (numRows > 0 && ptr[numRows] == 0) {
    // All the off-diagonal entries of A are zero, and all the
    // diagonal entries are (implicitly) 1.  Therefore compute: X :=
    // (1 - omega) X + omega B.  There's no need to care about the
    // direction or row ordering, since there are no cross-row data
    // dependencies in this case.
    const MatrixScalar oneMinusOmega = ArithTraits<MatrixScalar>::one() - omega;
    for (OffsetType j = 0; j < theNumCols; ++j) {
      RangeScalar* const x_j        = X + j * x_stride;
      const DomainScalar* const b_j = B + j * b_stride;
      for (OffsetType i = 0; i < theNumRows; ++i) {
        x_j[i] = oneMinusOmega * x_j[i] + omega * b_j[i];
      }
    }
    return;
  }

  if (numCols == 1) {
    if (direction[0] == 'F' || direction[0] == 'f') {
      for (LO ii = 0; ii < numRowInds; ++ii) {
        LO i               = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    } else if (direction[0] == 'B' || direction[0] == 'b') {
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO ii = numRowInds - 1; ii != 0; --ii) {
        LO i               = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
      {  // last loop iteration
        const LO ii        = 0;
        LO i               = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero();
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    }
  } else {  // numCols > 1
    // mfh 20 Dec 2012: If Gauss-Seidel for multivectors with
    // multiple columns becomes important, we can add unrolled
    // implementations.  The implementation below is not unrolled.
    // It may also be reasonable to parallelize over right-hand
    // sides, if there are enough of them, especially if the matrix
    // fits in cache.
    std::vector<RangeScalar> temp(numCols);
    RangeScalar* const x_temp = numCols == 0 ? NULL : &temp[0];

    if (direction[0] == 'F' || direction[0] == 'f') {
      for (LO ii = 0; ii < numRowInds; ++ii) {
        LO i = rowInd[ii];
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = Kokkos::ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
    } else if (direction[0] == 'B' || direction[0] == 'b') {  // backward mode
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO ii = numRowInds - 1; ii != 0; --ii) {
        LO i = rowInd[ii];
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = Kokkos::ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
      {  // last loop iteration
        const LO ii = 0;
        LO i        = rowInd[ii];
        for (OffsetType c = 0; c < theNumCols; ++c) {
          x_temp[c] = Kokkos::ArithTraits<RangeScalar>::zero();
        }
        for (OffsetType k = ptr[i]; k < ptr[i + 1]; ++k) {
          const LO j              = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < theNumCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride * c];
          }
        }
        for (OffsetType c = 0; c < theNumCols; ++c) {
          X[i + x_stride * c] += omega * D[i] * (B[i + b_stride * c] - x_temp[c]);
        }
      }
    }
  }
}

}  // namespace Sequential
}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_IMPL_SOR_HPP

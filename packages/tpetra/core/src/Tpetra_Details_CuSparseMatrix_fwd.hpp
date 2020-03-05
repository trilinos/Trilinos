#ifndef TPETRA_DETAILS_CUSPARSEMATRIX_FWD_HPP
#define TPETRA_DETAILS_CUSPARSEMATRIX_FWD_HPP

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_CUSPARSE

#include "Tpetra_Details_DefaultTypes.hpp"
#include <memory>

namespace Tpetra {
namespace Details {

/// \class CuSparseMatrix
/// \brief Opaque wrapper for sparse matrix arguments to cuSPARSE
///   functions.
///
/// \note To developers: Do not expose the declaration of this class
///   to downstream code.  Users should only deal with this class by
///   the forward declaration and functions available in this header
///   file.  Do not expose cuSPARSE headers or extern declarations to
///   downstream code.
class CuSparseMatrix;

namespace Impl {

/// \brief Delete the CuSparseMatrix.
///
/// This exists so that we can hide the declaration of CuSparseMatrix.
void deleteCuSparseMatrix(CuSparseMatrix* p);

} // namespace Impl

/// \fn getCuSparseMatrix
/// \brief Get cuSPARSE sparse matrix wrapper corresponding to the
///   given input data.
///
/// Input parameters correspond to exactly the parameters that
/// cusparseXcsrmv (old interface, deprecated in CUDA 10.1) and
/// cusparseCreateCsr (new interface, introduced in CUDA 10.1) both
/// expect.  This means that callers are responsible for making and
/// storing a copy of \c ptr with the correct element type, as would
/// be required if using this interface with
/// Tpetra::CrsMatrix::local_matrix_type.
///
/// \warning Do not let the returned object persist beyond
///   Kokkos::finalize.
///
/// \note We don't support complex types here, because Kokkos::complex
///   and CUDA's internal complex type do not have the same alignment
///   requirements.

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  float* val);

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  double* val);

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEMATRIX_FWD_HPP

#ifndef TPETRA_DETAILS_CUSPARSEVECTOR_FWD_HPP
#define TPETRA_DETAILS_CUSPARSEVECTOR_FWD_HPP

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_CUSPARSE

#include <cstdint> // int64_t
#include <memory>

namespace Kokkos {
class Cuda;
} // namespace Kokkos

namespace Tpetra {
namespace Details {

/// \class CuSparseVector
/// \brief Opaque wrapper for dense vector arguments to cuSPARSE
///   functions.
///
/// \note To developers: Do not expose the declaration of this class
///   to downstream code.  Users should only deal with this class by
///   the forward declaration and functions available in this header
///   file.  Do not expose cuSPARSE headers or extern declarations to
///   downstream code.
class CuSparseVector;

namespace Impl {

/// \brief Delete the CuSparseVector.
///
/// This exists so that we can hide the declaration of CuSparseVector.
void deleteCuSparseVector(CuSparseVector* p);

} // namespace Impl

/// \fn getCuSparseVector
/// \brief Get cuSPARSE vector wrapper corresponding to the given
///   input data.
///
/// \param data [in] Device pointer to the vector's data.
/// \param numRows [in] Number of local entries (rows).
///
/// \warning Do not let the returned object persist beyond
///   Kokkos::finalize.
///
/// \note We don't support complex types here, because Kokkos::complex
///   and CUDA's internal complex type do not have the same alignment
///   requirements.

std::unique_ptr<CuSparseVector, decltype(&Impl::deleteCuSparseVector)>
getCuSparseVector(float* const data, const int64_t numRows);

std::unique_ptr<CuSparseVector, decltype(&Impl::deleteCuSparseVector)>
getCuSparseVector(double* const data, const int64_t numRows);

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEVECTOR_FWD_HPP

#ifndef TPETRA_DETAILS_CUSPARSEHANDLE_FWD_HPP
#define TPETRA_DETAILS_CUSPARSEHANDLE_FWD_HPP

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_CUSPARSE

#include <memory>

namespace Kokkos {
class Cuda;
} // namespace Kokkos

namespace Tpetra {
namespace Details {

/// \class CuSparseHandle
/// \brief Opaque wrapper for cusparseHandle_t
///   (cuSPARSE handle instance)
///
/// \note To developers: Do not expose the declaration of this class
///   to downstream code.  Users should only deal with this class by
///   the forward declaration and functions available in this header
///   file.  Do not expose cuSPARSE headers or extern declarations to
///   downstream code.
class CuSparseHandle;

/// \brief Get cuSPARSE handle wrapper corresponding to the given
///   execution space instance (CUDA stream).
///
/// \param execSpace [in] Execution space instance (CUDA stream)
///   on which to run.
///
/// Tpetra may cache the result of the first call.  If the input's
/// CUDA stream is the same as that result's CUDA stream, then Tpetra
/// may return that cached result.
///
/// \warning Do not let the returned CuSparseHandle persist beyond
///   Kokkos::finalize.
///
/// \warning Kokkos::Cuda does not (currently) own its CUDA stream, so
///   you, the caller, are responsible for ensuring that the CUDA
///   stream lives at least as long as the returned CuSparseHandle.
///   This will always be the case if you use the default execution
///   space instance (<tt>Kokkos::Cuda()</tt>).
std::shared_ptr<CuSparseHandle>
getCuSparseHandle(const Kokkos::Cuda& execSpace);

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEHANDLE_FWD_HPP

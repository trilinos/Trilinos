#ifndef TPETRA_DETAILS_CUSPARSEHANDLE_HPP
#define TPETRA_DETAILS_CUSPARSEHANDLE_HPP

#include "Tpetra_Details_CuSparseHandle_fwd.hpp"
#ifdef HAVE_TPETRACORE_CUSPARSE
#include <cusparse.h>

// CUDA 10.1 introduced a new cuSPARSE interface and deprecated the
// old one.  CUDART_VERSION / 1000 gives the major release number (10,
// in this case); (CUDART_VERSION % 100) / 10 gives the minor release
// number (1, in this case).

#if defined(CUDART_VERSION) && CUDART_VERSION >= 10010
#  define HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
#endif

namespace Tpetra {
namespace Details {

class CuSparseHandle {
public:
  CuSparseHandle() = delete;
  CuSparseHandle(const CuSparseHandle&) = delete;
  CuSparseHandle& operator= (const CuSparseHandle&) = delete;
  CuSparseHandle(CuSparseHandle&&) = delete;
  CuSparseHandle& operator= (CuSparseHandle&&) = delete;

  CuSparseHandle(cusparseHandle_t handle,
                 const bool owning);
  ~CuSparseHandle();
  cusparseHandle_t getHandle() const;

private:
  // Kokkos::Cuda (the execution space) doesn't currently own its CUDA
  // stream, so we don't need to keep the instance given to the
  // constructor.  Also, cusparseHandle_t stores its CUDA stream; we
  // can get it later by calling cusparseGetStream.
  //
  // cusparseHandle_t is actually a pointer type.
  cusparseHandle_t handle_ {nullptr};
  bool owning_ {false};
};

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEHANDLE_HPP

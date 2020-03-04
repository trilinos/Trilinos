#ifndef TPETRA_DETAILS_CUSPARSEVECTOR_HPP
#define TPETRA_DETAILS_CUSPARSEVECTOR_HPP

#include "Tpetra_Details_CuSparseVector_fwd.hpp"
#ifdef HAVE_TPETRACORE_CUSPARSE
#include <cusparse.h>

// CUDA 10.1 introduced a new cuSPARSE interface and deprecated the
// old one.  CUDART_VERSION / 1000 gives the major release number (10,
// in this case); (CUDART_VERSION % 100) / 10 gives the minor release
// number (1, in this case).

#if defined(CUDART_VERSION) && CUDART_VERSION > 10010
#  define HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
#endif

namespace Tpetra {
namespace Details {

namespace Impl {

cudaDataType getCudaDataType(const float*);
cudaDataType getCudaDataType(const double*);

} // namespace Impl

// This class is always responsible for freeing its handle.
class CuSparseVector {
public:
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  // This is actually a pointer type.
  using handle_type = cusparseDnVecDescr_t;
#else
  struct vector_type {
    vector_type(int64_t size_, void* values_, cudaDataType valueType_) :
      size(size_), values(values_), valueType(valueType_)
    {}

    int64_t size = 0;
    void* values = nullptr;
    cudaDataType valueType = CUDA_R_64F;
  };
  using handle_type = vector_type*;
#endif

  CuSparseVector() = delete;
  CuSparseVector(const CuSparseVector&) = delete;
  CuSparseVector& operator= (const CuSparseVector&) = delete;
  CuSparseVector(CuSparseVector&&) = delete;
  CuSparseVector& operator= (CuSparseVector&&) = delete;

  CuSparseVector(int64_t size,
                 void* values,
                 cudaDataType valueType);
  ~CuSparseVector();

  handle_type getHandle() const;

private:
  handle_type handle_ {nullptr};
};

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEVECTOR_HPP

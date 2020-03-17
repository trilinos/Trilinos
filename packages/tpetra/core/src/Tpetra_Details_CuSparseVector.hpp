#ifndef TPETRA_DETAILS_CUSPARSEVECTOR_HPP
#define TPETRA_DETAILS_CUSPARSEVECTOR_HPP

#include "Tpetra_Details_CuSparseVector_fwd.hpp"
#ifdef HAVE_TPETRACORE_CUSPARSE
#include <cusparse.h>
// This include gives us HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE.
#include "Tpetra_Details_CuSparseHandle.hpp"

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

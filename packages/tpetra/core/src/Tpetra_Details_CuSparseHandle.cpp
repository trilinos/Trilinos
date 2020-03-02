#include "Tpetra_Details_CuSparseHandle.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include <cuda_runtime.h>

namespace Tpetra {
namespace Details {

cusparseHandle_t cuSparseRawHandle_ = nullptr;

namespace Impl {

void
createRawHandleAndSetStream(cusparseHandle_t& handle,
                            cudaStream_t stream)
{
  auto status = cusparseCreate(&handle);
  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  if (stream != nullptr) {
    status = cusparseSetStream(handle, stream);
    TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  }
}

std::shared_ptr<CuSparseHandle>
createCuSparseHandleSingleton(const Kokkos::Cuda& execSpace)
{
  auto finalizer = [] () {
    if (cuSparseRawHandle_ != nullptr) {
      (void) cusparseDestroy(cuSparseRawHandle_);
      cuSparseRawHandle_ = nullptr;
    }
  };
  Kokkos::push_finalize_hook(finalizer);
  createRawHandleAndSetStream(cuSparseRawHandle_,
                              execSpace.cuda_stream());
  // The Kokkos finalize hook will delete the handle.
  const bool owning = false;
  return std::shared_ptr<CuSparseHandle>
    (new CuSparseHandle(cuSparseRawHandle_, owning));
}

} // namespace Impl

CuSparseHandle::
CuSparseHandle(cusparseHandle_t handle, const bool owning) :
  handle_(handle), owning_(owning)
{}

CuSparseHandle::~CuSparseHandle()
{
  if (owning_ && handle_ != nullptr) {
    (void) cusparseDestroy(handle_);
  }
}

cusparseHandle_t
CuSparseHandle::getHandle() const {
  return handle_;
}

std::shared_ptr<CuSparseHandle>
getCuSparseHandle(const Kokkos::Cuda& execSpace)
{
  static std::shared_ptr<CuSparseHandle> singleton_;
  if (singleton_.get() == nullptr) {
    singleton_ = Impl::createCuSparseHandleSingleton(execSpace);
  }
  else {
    TEUCHOS_ASSERT( cuSparseRawHandle_ != nullptr );
    cudaStream_t curStrm = nullptr;
    auto status = cusparseGetStream(cuSparseRawHandle_, &curStrm);
    TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
    cudaStream_t inStrm = execSpace.cuda_stream();

    if (curStrm != inStrm) {
      // Input stream differs from singleton's stream, so don't use
      // the singleton; instead, create a new owning handle.
      cusparseHandle_t rawHandle = nullptr;
      Impl::createRawHandleAndSetStream(rawHandle, inStrm);
      const bool owning = true;
      return std::shared_ptr<CuSparseHandle>
        (new CuSparseHandle(rawHandle, owning));
    }
  }

  TEUCHOS_ASSERT( cuSparseRawHandle_ != nullptr );
  TEUCHOS_ASSERT( singleton_.get() != nullptr );
  return singleton_;
}

} // namespace Details
} // namespace Tpetra
#endif // HAVE_TPETRACORE_CUSPARSE

#include "Tpetra_Details_CuSparseVector.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

void deleteCuSparseVector(CuSparseVector* p)
{
  if (p != nullptr) {
    delete p;
  }
}

cudaDataType getCudaDataType(const float*) {
  return CUDA_R_32F;
}

cudaDataType getCudaDataType(const double*) {
  return CUDA_R_64F;
}
} // namespace Impl

CuSparseVector::
CuSparseVector(int64_t size,
               void* values,
               cudaDataType valueType)
{
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  const char prefix[] = "Tpetra::Details::CuSparseVector::"
    "CuSparseVector: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (values == nullptr && size != 0, std::invalid_argument,
     prefix << "size=" << size << " != 0, but values is null.");

  // With CUDA 10.1, cusparseCreateDnVec treats size=0 as an error,
  // even if values is nullptr.  Thus, it's best not to create the
  // handle in that case.
  if (size == 0) {
    handle_ = nullptr;
  }
  else {
    const cusparseStatus_t status =
      cusparseCreateDnVec(&handle_, size, values, valueType);
    TEUCHOS_TEST_FOR_EXCEPTION
      (status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
       prefix << "cusparseCreateDnVec(&handle_, size=" << size
       << ", values=" << values << ", valueType) returned an "
       "error code not CUSPARSE_STATUS_SUCCESS.");
  }
#else
  handle_ = new vector_type{size, values, valueType};
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
}

CuSparseVector::~CuSparseVector()
{
  if (handle_ != nullptr) {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    (void) cusparseDestroyDnVec(handle_);
#else
    delete handle_;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
}

CuSparseVector::handle_type
CuSparseVector::getHandle() const {
  return handle_;
}

namespace Impl {
template<class Scalar>
std::unique_ptr<CuSparseVector, decltype(&Impl::deleteCuSparseVector)>
getCuSparseVectorGeneric(Scalar* const data,
                         const int64_t numRows)
{
  const auto valueType = getCudaDataType(data);

  // make_unique is a C++14 feature.
  //
  // return std::make_unique<CuSparseVector>(numRows, data, valueType);

  return {new CuSparseVector(numRows, data, valueType),
          &Impl::deleteCuSparseVector};
}
} // namespace Impl

std::unique_ptr<CuSparseVector, decltype(&Impl::deleteCuSparseVector)>
getCuSparseVector(float* const data,
                  const int64_t numRows)
{
  return Impl::getCuSparseVectorGeneric<float>(data, numRows);
}

std::unique_ptr<CuSparseVector, decltype(&Impl::deleteCuSparseVector)>
getCuSparseVector(double* const data,
                  const int64_t numRows)
{
  return Impl::getCuSparseVectorGeneric<double>(data, numRows);
}

} // namespace Details
} // namespace Tpetra
#endif // HAVE_TPETRACORE_CUSPARSE

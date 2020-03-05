#include "Tpetra_Details_CuSparseMatrix.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Tpetra_Details_CuSparseHandle.hpp"
#include "Tpetra_Details_CuSparseVector.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>

namespace Tpetra {
namespace Details {
namespace Impl {

#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
cusparseIndexType_t getIndexType(const int32_t*) {
  return CUSPARSE_INDEX_32I;
}
cusparseIndexType_t getIndexType(const int64_t*) {
  return CUSPARSE_INDEX_64I;
}
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

void deleteCuSparseMatrix(CuSparseMatrix* p)
{
  if (p != nullptr) {
    delete p;
  }
}

} // namespace Impl

CuSparseMatrix::
CuSparseMatrix(const int64_t numRows,
               const int64_t numCols,
               const int64_t numEntries,
               local_ordinal_type* ptr,
               local_ordinal_type* ind,
               void* val,
               const cudaDataType valueType)
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  : valueType_(valueType)
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
{
  cusparseStatus_t status;

#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  const cusparseIndexType_t csrRowOffsetsType =
    Impl::getIndexType(ptr);
  const cusparseIndexType_t csrColIndType =
    Impl::getIndexType(ind);
  status =
    cusparseCreateCsr(&handle_, numRows, numCols, numEntries,
                      ptr, ind, val,
                      csrRowOffsetsType, csrColIndType,
                      idxBase, valueType);
  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
#else
  handle_ = new matrix_type{numRows, numCols, numEntries,
                            ptr, ind, val, valueType};

  status = cusparseCreateMatDescr(&descr_);
  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );

  try {
    status = cusparseSetMatType(descr_, CUSPARSE_MATRIX_TYPE_GENERAL);
    TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  }
  catch (...) {
    freeHandle(handle_);
    throw;
  }

  try {
    // GENERAL matrix type doesn't use fill mode
    // status = cusparseSetMatFillMode(descr_, whatever);

    status = cusparseSetMatDiagType(descr_, CUSPARSE_DIAG_TYPE_NON_UNIT);
    TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );

    status = cusparseSetMatIndexBase(descr_, CUSPARSE_INDEX_BASE_ZERO);
    TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  }
  catch (...) {
    freeHandle(handle_);
    freeDescr(descr_);
    throw;
  }
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
}

void
CuSparseMatrix::freeBuffer(void* buffer)
{
  if (buffer != nullptr) {
    (void) cudaFree(buffer);
  }
}

void
CuSparseMatrix::freeHandle(handle_type handle)
{
  if (handle != nullptr) {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    (void) cusparseDestroySpMat(handle);
#else
    delete handle;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
}

#ifndef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
void
CuSparseMatrix::freeDescr(cusparseMatDescr_t descr)
{
  if (descr != nullptr) {
    (void) cusparseDestroyMatDescr(descr);
  }
}

cusparseMatDescr_t
CuSparseMatrix::getDescr() const {
  return descr_;
}
#endif // NOT HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

CuSparseMatrix::~CuSparseMatrix()
{
  freeHandle(handle_);
#ifndef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  freeDescr(descr_);
#endif // NOT HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  freeBuffer(buffer_);
}

CuSparseMatrix::handle_type
CuSparseMatrix::getHandle() const {
  return handle_;
}

cudaDataType
CuSparseMatrix::getValueType() const {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  return valueType_;
#else
  return handle_->valueType;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
}

CuSparseMatrix::algorithm_t
CuSparseMatrix::getAlgorithm(const cusparseOperation_t op) const {
  if (op == CUSPARSE_OPERATION_NON_TRANSPOSE) {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    return CUSPARSE_CSRMV_ALG2; // merge path
#else
    return CUSPARSE_ALG_MERGE_PATH;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
  else {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    return CUSPARSE_MV_ALG_DEFAULT;
#else
    return CUSPARSE_ALG_NAIVE;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
}

void*
CuSparseMatrix::
reallocBufferIfNeededAndGetBuffer(const size_t minNeededBufSize)
{
  // FIXME (mfh 05 Mar 2020) Hack, since cusparseCsrmvEx_bufferSize
  // claims a buffer size of -1 even though its status is
  // CUSPARSE_STATUS_SUCCESS.  cuSPARSE requires the buffer to have
  // correct alignment.  We're only calling cuSPARSE for float and
  // double, but it's no more expensive to allocate a tiny bit more
  // for cuComplexDouble alignment.
  if (minNeededBufSize == static_cast<size_t>(-1)) {
    return reallocBufferIfNeededAndGetBuffer(16);
  }

  if (minNeededBufSize > bufSize_) {
    freeBuffer(buffer_);
    buffer_ = nullptr;
    const cudaError_t status = cudaMalloc(&buffer_, minNeededBufSize);

    if (status != cudaSuccess) {
      const std::string statusStr = status == cudaErrorInvalidValue ?
        "cudaErrorInvalidValue" :
        (status == cudaErrorMemoryAllocation ?
         "cudaErrorMemoryAllocation" :
         "UNKNOWN");
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::runtime_error, "cudaMalloc(&buffer_, "
         << minNeededBufSize << ") returned status=" << statusStr
         << ".");
    }
  }
  return buffer_;
}

namespace Impl {
template<class Scalar>
std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrixGeneric(const int64_t numRows,
                         const int64_t numCols,
                         const int64_t numEntries,
                         DefaultTypes::local_ordinal_type* ptr,
                         DefaultTypes::local_ordinal_type* ind,
                         Scalar* val)
{
  const auto valueType = getCudaDataType(val);
  return {new CuSparseMatrix(numRows, numCols, numEntries,
                             ptr, ind, val, valueType),
          &Impl::deleteCuSparseMatrix};
}

cusparseOperation_t
getCuSparseOperation(const Teuchos::ETransp mode)
{
  if (mode == Teuchos::NO_TRANS) {
    return CUSPARSE_OPERATION_NON_TRANSPOSE;
  }
  else if (mode == Teuchos::TRANS) {
    return CUSPARSE_OPERATION_TRANSPOSE;
  }
  else { // mode == Teuchos::CONJ_TRANS
    return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
  }
}

template<class Scalar>
void
cuSparseMatrixVectorMultiplyGeneric(
  CuSparseHandle& handle,
  const Teuchos::ETransp operation,
  const Scalar alpha,
  CuSparseMatrix& matrix,
  CuSparseVector& x,
  const Scalar beta,
  CuSparseVector& y)
{
  cusparseHandle_t rawHandle = handle.getHandle();
  cusparseOperation_t rawOp = getCuSparseOperation(operation);
  auto rawMatrix = matrix.getHandle();
  cudaDataType valType = matrix.getValueType();
  auto alg = matrix.getAlgorithm(rawOp);

  size_t minNeededBufSize = 0;
  cusparseStatus_t status;
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  status =
    cusparseSpMV_bufferSize(rawHandle, rawOp, &alpha, rawMatrix,
                            x.getHandle(), &beta, y.getHandle(),
                            valType, alg, &minNeededBufSize);
  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  void* buf =
    matrix.reallocBufferIfNeededAndGetBuffer(minNeededBufSize);
  status =
    cusparseSpMV(rawHandle, rawOp, &alpha, rawMatrix, x.getHandle(),
                 &beta, y.getHandle(), valType, alg, buf);
#else
  cusparseMatDescr_t descrA = matrix.getDescr();
  const int m = rawMatrix->numRows;
  const int n = rawMatrix->numCols;
  const int nnz = rawMatrix->numEntries;
  const Scalar* val = reinterpret_cast<const Scalar*>(rawMatrix->val);
  const int* ptr = reinterpret_cast<const int*>(rawMatrix->ptr);
  const int* ind = reinterpret_cast<const int*>(rawMatrix->ind);

  // This is the "data type used for computation."  CUDA 10.1 doesn't
  // currently support this being different than val's or x's
  // cudaDataType, if the latter are CUDA_R_32F or CUDA_R_64F.  That's
  // a bit unfortunate for float, but ah well.  If that ever should
  // change, consider making execType=CUDA_R_64F when Scalar=float.
  cudaDataType execType = valType;

  status = cusparseCsrmvEx_bufferSize
    (rawHandle, matrix.getAlgorithm(rawOp), rawOp, m, n, nnz,
     &alpha, valType, descrA, val, valType, ptr, ind,
     x.getHandle()->values, valType, &beta, valType,
     y.getHandle()->values, valType, execType, &minNeededBufSize);
  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
  void* buf =
    matrix.reallocBufferIfNeededAndGetBuffer(minNeededBufSize);
  status =
    cusparseCsrmvEx(rawHandle, matrix.getAlgorithm(rawOp), rawOp,
                    m, n, nnz, &alpha, valType, descrA, val, valType,
                    ptr, ind, x.getHandle()->values, valType,
                    &beta, valType, y.getHandle()->values, valType,
                    execType, buf);
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

  TEUCHOS_ASSERT( status == CUSPARSE_STATUS_SUCCESS );
}

} // namespace Impl

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  float* val)
{
  return Impl::getCuSparseMatrixGeneric<float>
    (numRows, numCols, numEntries, ptr, ind, val);
}

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  double* val)
{
  return Impl::getCuSparseMatrixGeneric<double>
    (numRows, numCols, numEntries, ptr, ind, val);
}

void
cuSparseMatrixVectorMultiply(
  CuSparseHandle& handle,
  const Teuchos::ETransp operation,
  const double alpha,
  CuSparseMatrix& matrix,
  CuSparseVector& x,
  const double beta,
  CuSparseVector& y)
{
  Impl::cuSparseMatrixVectorMultiplyGeneric<double>
    (handle, operation, alpha, matrix, x, beta, y);
}

void
cuSparseMatrixVectorMultiply(
  CuSparseHandle& handle,
  const Teuchos::ETransp operation,
  const float alpha,
  CuSparseMatrix& matrix,
  CuSparseVector& x,
  const float beta,
  CuSparseVector& y)
{
  Impl::cuSparseMatrixVectorMultiplyGeneric<double>
    (handle, operation, alpha, matrix, x, beta, y);
}

} // namespace Details
} // namespace Tpetra
#endif // HAVE_TPETRACORE_CUSPARSE

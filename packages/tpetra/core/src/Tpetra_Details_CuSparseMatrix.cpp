#include "Tpetra_Details_CuSparseMatrix.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_CuSparseHandle.hpp"
#include "Tpetra_Details_CuSparseVector.hpp"
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include <sstream>
#include <stdexcept>

#define TPETRA_DETAILS_CUSPARSE_CHECK_STATUS( status ) \
  TEUCHOS_TEST_FOR_EXCEPTION \
    (status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, \
     "cuSPARSE call returned status = " \
     << Impl::cuSparseStatusString(status) << ".");

namespace Tpetra {
namespace Details {
namespace Impl {

std::string cuSparseStatusString(cusparseStatus_t status)
{
  if (status == CUSPARSE_STATUS_SUCCESS) {
    return "CUSPARSE_STATUS_SUCCESS";
  }
  else if (status == CUSPARSE_STATUS_NOT_INITIALIZED) {
    return "CUSPARSE_STATUS_NOT_INITIALIZED";
  }
  else if (status == CUSPARSE_STATUS_ALLOC_FAILED) {
    return "CUSPARSE_STATUS_ALLOC_FAILED";
  }
  else if (status == CUSPARSE_STATUS_INVALID_VALUE) {
    return "CUSPARSE_STATUS_INVALID_VALUE";
  }
  else if (status == CUSPARSE_STATUS_ARCH_MISMATCH) {
    return "CUSPARSE_STATUS_ARCH_MISMATCH";
  }
  else if (status == CUSPARSE_STATUS_EXECUTION_FAILED) {
    return "CUSPARSE_STATUS_EXECUTION_FAILED";
  }
  else if (status == CUSPARSE_STATUS_INTERNAL_ERROR) {
    return "CUSPARSE_STATUS_INTERNAL_ERROR";
  }
  else if (status == CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED) {
    return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
  }
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  else if (status == CUSPARSE_STATUS_NOT_SUPPORTED) {
    return "CUSPARSE_STATUS_MAPPING_ERROR";
  }
#else
  else if (status == CUSPARSE_STATUS_MAPPING_ERROR) {
    return "CUSPARSE_STATUS_MAPPING_ERROR";
  }
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  else {
    return "UNKNOWN (this is bad)";
  }
}

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
               const cudaDataType valueType,
               const CuSparseMatrixVectorMultiplyAlgorithm alg) :
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  valueType_(valueType),
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  alg_(alg)
{
  using std::endl;
  const char funcName[] = "Tpetra::Details::CuSparseMatrix::"
    "CuSparseMatrix";

  const bool debug = Details::Behavior::debug("cuSPARSE");
  const bool verbose = Details::Behavior::verbose("cuSPARSE");
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    std::ostringstream os;
    os << funcName << ": ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start" << endl;
    std::cerr << os.str();
  }
  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On entry to "
       "Tpetra::Details::CuSparseMatrix's constructor, CUDA is in an "
       "erroneous state \"" << cudaGetErrorName(lastErr) << "\".");
  }

  cusparseStatus_t status;
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  if (numEntries == 0) {
    // mfh 17 Mar 2020: numRows=0 or numCols=0 are fine, but CUDA 10.1
    // treats numEntries=0 as an error.
    handle_ = nullptr;
  }
  else {
    const cusparseIndexType_t csrRowOffsetsType =
      Impl::getIndexType(ptr);
    const cusparseIndexType_t csrColIndType =
      Impl::getIndexType(ind);
    const cusparseIndexBase_t idxBase = CUSPARSE_INDEX_BASE_ZERO;
    status =
      cusparseCreateCsr(&handle_, numRows, numCols, numEntries,
                        ptr, ind, val,
                        csrRowOffsetsType, csrColIndType,
                        idxBase, valueType);
    TEUCHOS_TEST_FOR_EXCEPTION
      (status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
       funcName << ": cusparseCreateCsr(&handle_, numRows=" << numRows
       << ", numCols=" << numCols << ", numEntries=" << numEntries
       << ", ptr=" << ptr << ", ind=" << ind << ", val=" << val
       << "...) returned status = "
       << Impl::cuSparseStatusString(status) << ".");
  }
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

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str();
  }
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
  if (alg_ == CuSparseMatrixVectorMultiplyAlgorithm::DEFAULT ||
      op != CUSPARSE_OPERATION_NON_TRANSPOSE) {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    return CUSPARSE_MV_ALG_DEFAULT;
#else
    return CUSPARSE_ALG_NAIVE;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
  else {
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
    return CUSPARSE_CSRMV_ALG2; // merge path
#else
    return CUSPARSE_ALG_MERGE_PATH;
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  }
}

#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
void*
CuSparseMatrix::
reallocBufferIfNeededAndGetBuffer(const size_t minNeededBufSize)
{
  // FIXME (mfh 05 Mar 2020) Hack, since cusparseCsrmvEx_bufferSize
  // might claim a buffer size of -1 even though its status is
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
#else
void*
CuSparseMatrix::
reallocBufferIfNeededAndGetBuffer(const size_t minNeededBufSize)
{
  using std::endl;
  const bool debug = Details::Behavior::debug("cuSPARSE");
  const bool verbose = Details::Behavior::verbose("cuSPARSE");
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    std::ostringstream os;
    os << "Tpetra::Details::CuSparseMatrix::"
      "reallocBufferIfNeededAndGetBuffer: ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start: minNeededBufSize=" << minNeededBufSize << endl;
    std::cerr << os.str();
  }
  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On entry to "
       "Tpetra::Details::CuSparseMatrix::"
       "reallocBufferIfNeededAndGetBuffer, CUDA is in an "
       "erroneous state \"" << cudaGetErrorName(lastErr) << "\".");
  }

  // NOTE (mfh 09 Mar 2020) cusparseCsrmvEx_bufferSize doesn't
  // actually use the buffer, but cuSPARSE's documentation claims that
  // cusparseCsrmvEx needs the buffer to have correct alignment for
  // the matrix's value type.  We're only calling cuSPARSE for float
  // and double, but it's no more expensive to allocate a tiny bit
  // more for cuComplexDouble alignment.
  constexpr size_t requiredBufSize = 16;

  if (buffer_ == nullptr) {
    const cudaError_t status = cudaMalloc(&buffer_, requiredBufSize);
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
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

namespace Impl {
template<class Scalar>
std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrixGeneric(const int64_t numRows,
                         const int64_t numCols,
                         const int64_t numEntries,
                         DefaultTypes::local_ordinal_type* ptr,
                         DefaultTypes::local_ordinal_type* ind,
                         Scalar* val,
                         const CuSparseMatrixVectorMultiplyAlgorithm alg)
{
  const auto valueType = getCudaDataType(val);
  return {new CuSparseMatrix(numRows, numCols, numEntries,
                             ptr, ind, val, valueType, alg),
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
  using std::endl;
  const char funcName[] =
    "Tpetra::Details::cuSparseMatrixVectorMultiplyGeneric";

  const bool debug = Details::Behavior::debug("cuSPARSE");
  const bool verbose = Details::Behavior::verbose("cuSPARSE");
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    std::ostringstream os;
    os << funcName << ": ";
    prefix = std::unique_ptr<std::string>(new std::string(os.str()));
    os << "Start" << endl;
    std::cerr << os.str();
  }
  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, funcName << ": On "
       "entry to this function, CUDA is in an erroneous state \"" <<
       cudaGetErrorName(lastErr) << "\".");
  }

  auto rawMatrix = matrix.getHandle();
  // mfh 17 Mar 2020: We should have already handled the case where
  // the matrix handle is null.  cuSPARSE's CUDA >= 10.1 interface
  // treats creating a matrix handle with zero entries as an error, so
  // we just make the matrix handle null in that case.  Similarly, the
  // CUDA >= 10.1 interface forbids creating dense vector handles for
  // vectors with zero entries, so we also make the vector handle null
  // in that case; we should have already handled that.
  TEUCHOS_ASSERT( rawMatrix != nullptr );

  cusparseHandle_t rawHandle = handle.getHandle();
  cusparseOperation_t rawOp = getCuSparseOperation(operation);
  cudaDataType valType = matrix.getValueType();
  auto alg = matrix.getAlgorithm(rawOp);

  size_t minNeededBufSize = 0;
  cusparseStatus_t status;
#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Call cusparseSpMV_bufferSize" << endl;
    std::cerr << os.str();
  }
  status =
    cusparseSpMV_bufferSize(rawHandle, rawOp, &alpha, rawMatrix,
                            x.getHandle(), &beta, y.getHandle(),
                            valType, alg, &minNeededBufSize);
  TEUCHOS_TEST_FOR_EXCEPTION
    (status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
     funcName << ": cusparseSpMV_bufferSize returned status = "
     << Impl::cuSparseStatusString(status) << ".");
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Buffer size: " << minNeededBufSize << endl;
    std::cerr << os.str();
  }
  void* buf =
    matrix.reallocBufferIfNeededAndGetBuffer(minNeededBufSize);
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Call cusparseSpMV" << endl;
    std::cerr << os.str();
  }
  status =
    cusparseSpMV(rawHandle, rawOp, &alpha, rawMatrix, x.getHandle(),
                 &beta, y.getHandle(), valType, alg, buf);
  TEUCHOS_TEST_FOR_EXCEPTION
    (status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
     funcName << ": cusparseSpMV returned status = "
     << Impl::cuSparseStatusString(status) << ".");

#else

  cusparseMatDescr_t descrA = matrix.getDescr();
  const int m = rawMatrix->numRows;
  const int n = rawMatrix->numCols;
  const int nnz = rawMatrix->numEntries;
  const Scalar* val = reinterpret_cast<const Scalar*>(rawMatrix->val);
  const int* ptr = reinterpret_cast<const int*>(rawMatrix->ptr);
  const int* ind = reinterpret_cast<const int*>(rawMatrix->ind);

  if (verbose) {
    cusparsePointerMode_t ptrMode;
    status = cusparseGetPointerMode(rawHandle, &ptrMode);
    TPETRA_DETAILS_CUSPARSE_CHECK_STATUS( status );
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "cuSPARSE pointer mode: ";
      if (ptrMode == CUSPARSE_POINTER_MODE_HOST) {
        os << "CUSPARSE_POINTER_MODE_HOST";
      }
      else if (ptrMode == CUSPARSE_POINTER_MODE_DEVICE) {
        os << "CUSPARSE_POINTER_MODE_DEVICE";
      }
      else {
        os << "UNKNOWN (this is bad)";
      }
      os << endl << *prefix << "m=" << m << ", n=" << n
         << ", nnz=" << nnz << ", x.getHandle()->size="
         << x.getHandle()->size << ", y.getHandle()->size="
         << y.getHandle()->size << endl;
      std::cerr << os.str();
    }
  }

  TEUCHOS_ASSERT( descrA != nullptr );
  TEUCHOS_ASSERT( rawHandle != nullptr );
  if (m != 0 && n != 0 && nnz != 0) {
    TEUCHOS_ASSERT( val != nullptr );
    TEUCHOS_ASSERT( ptr != nullptr );
    TEUCHOS_ASSERT( ind != nullptr );
  }
  if (m != 0 && nnz != 0) {
    TEUCHOS_ASSERT( y.getHandle()->values != nullptr );
  }
  if (n != 0 && nnz != 0) {
    TEUCHOS_ASSERT( x.getHandle()->values != nullptr );
  }

  // This is the "data type used for computation."  CUDA 10.1 doesn't
  // currently support this being different than val's or x's
  // cudaDataType, if the latter are CUDA_R_32F or CUDA_R_64F.  That's
  // a bit unfortunate for float, but ah well.  If that ever should
  // change, consider making execType=CUDA_R_64F when Scalar=float.
  cudaDataType execType = valType;

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Call cusparseCsrmvEx_bufferSize" << endl;
    std::cerr << os.str();
  }
  status = cusparseCsrmvEx_bufferSize
    (rawHandle, alg, rawOp, m, n, nnz,
     &alpha, valType, descrA, val, valType, ptr, ind,
     x.getHandle()->values, valType, &beta, valType,
     y.getHandle()->values, valType, execType, &minNeededBufSize);
  TPETRA_DETAILS_CUSPARSE_CHECK_STATUS( status );
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Buffer size: " << minNeededBufSize << endl;
    std::cerr << os.str();
  }
  void* buf =
    matrix.reallocBufferIfNeededAndGetBuffer(minNeededBufSize);
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Call cusparseCsrmvEx" << endl;
    std::cerr << os.str();
  }
  status =
    cusparseCsrmvEx(rawHandle, alg, rawOp,
                    m, n, nnz, &alpha, valType, descrA, val, valType,
                    ptr, ind, x.getHandle()->values, valType,
                    &beta, valType, y.getHandle()->values, valType,
                    execType, buf);
  TPETRA_DETAILS_CUSPARSE_CHECK_STATUS( status );
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

  if (debug) {
    const cudaError_t lastErr = cudaGetLastError();
    TEUCHOS_TEST_FOR_EXCEPTION
      (lastErr != cudaSuccess, std::runtime_error, "On entry to "
       "Tpetra::Details::cuSparseMatrixVectorMultiplyGeneric, CUDA "
       "is in an erroneous state \"" << cudaGetErrorName(lastErr)
       << "\".");
  }
  if (verbose) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str();
  }
}

} // namespace Impl

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  float* val,
                  const CuSparseMatrixVectorMultiplyAlgorithm alg)
{
  return Impl::getCuSparseMatrixGeneric<float>
    (numRows, numCols, numEntries, ptr, ind, val, alg);
}

std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)>
getCuSparseMatrix(const int64_t numRows,
                  const int64_t numCols,
                  const int64_t numEntries,
                  DefaultTypes::local_ordinal_type* ptr,
                  DefaultTypes::local_ordinal_type* ind,
                  double* val,
                  const CuSparseMatrixVectorMultiplyAlgorithm alg)
{
  return Impl::getCuSparseMatrixGeneric<double>
    (numRows, numCols, numEntries, ptr, ind, val, alg);
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

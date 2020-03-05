#include "Tpetra_Details_CuSparseMatrix.hpp"

#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"

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
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

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

void
CuSparseMatrix::freeDescr(cusparseMatDescr_t descr)
{
  if (descr != nullptr) {
    (void) cusparseDestroyMatDescr(descr);
  }
}

CuSparseMatrix::~CuSparseMatrix()
{
  freeHandle(handle_);
  freeDescr(descr_);
}

CuSparseMatrix::handle_type
CuSparseMatrix::getHandle() const {
  return handle_;
}

cusparseMatDescr_t
CuSparseMatrix::getDescr() const {
  return descr_;
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

} // namespace Details
} // namespace Tpetra
#endif // HAVE_TPETRACORE_CUSPARSE

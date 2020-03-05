#ifndef TPETRA_DETAILS_CUSPARSEMATRIX_HPP
#define TPETRA_DETAILS_CUSPARSEMATRIX_HPP

#include "Tpetra_Details_CuSparseMatrix_fwd.hpp"
#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Tpetra_Details_CuSparseVector.hpp"

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
namespace Impl {

cusparseIndexType_t getIndexType(const int32_t*);
cusparseIndexType_t getIndexType(const int64_t*);

} // namespace Impl
#endif // HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE

// This class is always responsible for freeing its handle.
class CuSparseMatrix {
public:
  using local_ordinal_type =
    ::Tpetra::Details::DefaultTypes::local_ordinal_type;

#ifdef HAVE_TPETRACORE_CUSPARSE_NEW_INTERFACE
  // This is actually a pointer type.
  using handle_type = cusparseSpMatDescr_t;
#else
  struct matrix_type {
    matrix_type(const int64_t numRows_,
                const int64_t numCols_,
                const int64_t numEntries_,
                local_ordinal_type* ptr_,
                local_ordinal_type* ind_,
                void* val_,
                cudaDataType valueType_) :
      numRows(numRows_),
      numCols(numCols_),
      numEntries(numEntries_),
      ptr(ptr_),
      ind(ind_),
      val(val_),
      valueType(valueType_)
    {}

    int64_t numRows = 0;
    int64_t numCols = 0;
    int64_t numEntries = 0;
    local_ordinal_type* ptr = nullptr;
    local_ordinal_type* ind = nullptr;
    void* val = nullptr;
    // cusparseIndexType_t csrRowOffsetsType = Impl::getIndexType(ptr);
    // cusparseIndexType_t csrColIndType = Impl::getIndexType(ind);
    // cusparseIndexType_t idxBase = CUSPARSE_INDEX_BASE_ZERO;
    cudaDataType valueType = CUDA_R_64F;
  };
  using handle_type = matrix_type*;
#endif

  CuSparseMatrix() = delete;
  CuSparseMatrix(const CuSparseMatrix&) = delete;
  CuSparseMatrix& operator= (const CuSparseMatrix&) = delete;
  CuSparseMatrix(CuSparseMatrix&&) = delete;
  CuSparseMatrix& operator= (CuSparseMatrix&&) = delete;

  CuSparseMatrix(const int64_t numRows,
                 const int64_t numCols,
                 const int64_t numEntries,
                 local_ordinal_type* ptr,
                 local_ordinal_type* ind,
                 void* val,
                 const cudaDataType valueType);
  ~CuSparseMatrix();

  handle_type getHandle() const;

private:
  handle_type handle_ {nullptr};
};

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_CUSPARSE

#endif // TPETRA_DETAILS_CUSPARSEMATRIX_HPP

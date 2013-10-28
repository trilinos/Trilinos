//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

/// \file Kokkos_CUSPARSEOps.hpp
/// \brief CUSPARSE implementation of local sparse matrix kernels.
///
/// This header contains the CUSPARSE implementation of local sparse
/// matrix kernels. (the \c LocalMatOps fifth template parameter of
/// Tpetra::CrsMatrix).  "Local" means "within an MPI process."  This
/// is a KokkosClassic header file.

#ifndef KOKKOS_CUSPARSEOPS_HPP
#define KOKKOS_CUSPARSEOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Describable.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_CUDANodeUtils.hpp>

#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <cusparse_v2.h>

namespace Teuchos {
  template<>
  class TypeNameTraits<cusparseHandle_t> {
  public:
    static std::string name() { return std::string("cusparseHandle_t"); }
    static std::string concreteName( const cusparseHandle_t& ) {
      return std::string("cusparseHandle_t");
    }
  };
  template<>
  class TypeNameTraits<cusparseMatDescr_t> {
  public:
    static std::string name() { return std::string("cusparseMatDescr_t"); }
    static std::string concreteName( const cusparseMatDescr_t& ) {
      return std::string("cusparseMatDescr_t");
    }
  };
  template<>
  class TypeNameTraits<cusparseSolveAnalysisInfo_t> {
  public:
    static std::string name() { return std::string("cusparseSolveAnalysisInfo_t"); }
    static std::string concreteName( const cusparseSolveAnalysisInfo_t& ) {
      return std::string("cusparseSolveAnalysisInfo_t");
    }
  };
}

namespace KokkosClassic {

  namespace CUSPARSEdetails {

    class Session {
      public:
        //! Must be called to initialize the CUSPARSE session.
        static void init();

        //! Accessor for CUSPARSE session.
        static RCP<const cusparseHandle_t> getHandle();

      private:
        static RCP<cusparseHandle_t> session_handle_;
    };

    RCP<cusparseMatDescr_t>           createMatDescr();
    RCP<cusparseSolveAnalysisInfo_t>  createSolveAnalysisInfo();

    class CUSPARSESessionDestroyer {
    public:
      CUSPARSESessionDestroyer();
      void free(cusparseHandle_t *ptr);
    };

    class CUSPARSEMatDescDestroyer {
    public:
      CUSPARSEMatDescDestroyer();
      void free(cusparseMatDescr_t *ptr);
    };

    class CUSPARSESolveAnalysisDestroyer {
    public:
      CUSPARSESolveAnalysisDestroyer();
      void free(cusparseSolveAnalysisInfo_t *ptr);
    };


    template <typename T>
    struct CUSPARSEUnsupportedScalar
    {
      public:
      static inline cusparseStatus_t notSupported() {
        T::this_type_is_not_supported_by_CUSPARSE();
        return CUSPARSE_STATUS_ARCH_MISMATCH;
      }
    };


    template <class T>
    class CUSPARSETemplateAdaptors {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const T *alpha,
            const cusparseMatDescr_t descrA, const T *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const T *B, int ldb,
            const T *beta, T *C, int ldc)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const T *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const T *alpha,
                  const cusparseMatDescr_t descrA,
                  const T *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const T *X, int ldx,
                  T *Y, int ldy)
      { return CUSPARSEUnsupportedScalar<T>::notSupported(); }
    };

    // float support
#ifdef HAVE_KOKKOSCLASSIC_CUDA_FLOAT
    template <>
    class CUSPARSETemplateAdaptors<float>
    {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const float *alpha,
            const cusparseMatDescr_t descrA, const float *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const float *B, int ldb,
            const float *beta, float *C, int ldc)
      { return cusparseScsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const float *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseScsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const float *alpha,
                  const cusparseMatDescr_t descrA,
                  const float *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const float *X, int ldx,
                  float *Y, int ldy)
      { return cusparseScsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
    // double support
    template <>
    class CUSPARSETemplateAdaptors<double>
    {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const double *alpha,
            const cusparseMatDescr_t descrA, const double *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const double *B, int ldb,
            const double *beta, double *C, int ldc)
      { return cusparseDcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const double *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseDcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const double *alpha,
                  const cusparseMatDescr_t descrA,
                  const double *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const double *X, int ldx,
                  double *Y, int ldy)
      { return cusparseDcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT
    // complex<float> support
    template <>
    class CUSPARSETemplateAdaptors<std::complex<float> >
    {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const cuComplex *alpha,
            const cusparseMatDescr_t descrA, const cuComplex *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const cuComplex *B, int ldb,
            const cuComplex *beta, cuComplex *C, int ldc)
      { return cusparseCcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const cuComplex *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseCcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const cuComplex *alpha,
                  const cusparseMatDescr_t descrA,
                  const cuComplex *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const cuComplex *X, int ldx,
                  cuComplex *Y, int ldy)
      { return cusparseCcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE
    // complex<double> support
    template <>
    class CUSPARSETemplateAdaptors<std::complex<double> >
    {
      public:
      //
      static inline cusparseStatus_t
      CSRMM(cusparseHandle_t handle, cusparseOperation_t transA,
            int m, int n, int k, int nnz, const cuDoubleComplex *alpha,
            const cusparseMatDescr_t descrA, const cuDoubleComplex *csrValA,
            const int *csrRowPtrA, const int *csrColIndA,
            const cuDoubleComplex *B, int ldb,
            const cuDoubleComplex *beta, cuDoubleComplex *C, int ldc)
      { return cusparseZcsrmm(handle,transA,m,n,k,nnz,
                              alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                              B,ldb,beta,C,ldc); }
      //
      static inline cusparseStatus_t
      CSRSM_analysis(cusparseHandle_t handle, cusparseOperation_t transA,
                     int m, int nnz, const cusparseMatDescr_t descrA,
                     const cuDoubleComplex *csrValA, const int *csrRowPtrA,
                     const int *csrColIndA, cusparseSolveAnalysisInfo_t info)
      { return cusparseZcsrsm_analysis(handle,transA,m,nnz,descrA,
                                       csrValA,csrRowPtrA,csrColIndA,info); }
      //
      static inline cusparseStatus_t
      CSRSM_solve(cusparseHandle_t handle, cusparseOperation_t transA,
                  int m, int n, const cuDoubleComplex *alpha,
                  const cusparseMatDescr_t descrA,
                  const cuDoubleComplex *csrValA, const int *csrRowPtrA,
                  const int *csrColIndA, cusparseSolveAnalysisInfo_t info,
                  const cuDoubleComplex *X, int ldx,
                  cuDoubleComplex *Y, int ldy)
      { return cusparseZcsrsm_solve(handle,transA,m,n,
                                    alpha,descrA,csrValA,csrRowPtrA,csrColIndA,
                                    info,X,ldx,Y,ldy); }
    };
#endif

  } // end of namespace CUSPARSEdetails

  /// \class CUSPARSECrsGraph
  /// \brief Sparse graph class supporting the CUSPARSE library.
  /// \tparam Node The Kokkos Node type
  ///
  /// CUSPARSEOps uses this class to represent the local sparse graph.
  ///
  /// \note This class is really an implementation detail of Tpetra.
  ///   Tpetra users do not normally interact with this class.  In
  ///   fact, neither Tpetra::CrsGraph nor Tpetra::CrsMatrix give
  ///   users a way to access an instance of this class.
  template <class Node>
  class CUSPARSECrsGraph : public CrsGraphBase<int,Node> {
  public:
    CUSPARSECrsGraph (int numRows, int numCols, const RCP<Node> &node,
                      const RCP<ParameterList> &params);
    bool isEmpty() const;
    void setStructure(const ArrayRCP<const size_t>  &ptrs,
                      const ArrayRCP<const int> &inds);
    void setDeviceData(const ArrayRCP<const int> &devptrs,
                       const ArrayRCP<const int> &devinds);
    inline ArrayRCP<const size_t> getPointers() const;
    inline ArrayRCP<const int> getIndices() const;
    inline ArrayRCP<const int> getDevPointers() const;
    inline ArrayRCP<const int> getDevIndices() const;
    inline bool isInitialized() const;
    void setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag);
    void getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;
    RCP<cusparseMatDescr_t> getMatDesc() const;
  private:
    bool isInitialized_;
    bool isEmpty_;
    // cusparse matrix description handle
    RCP<cusparseMatDescr_t> matdescr_;
    Teuchos::EUplo uplo_;
    Teuchos::EDiag diag_;
    // graph data
    ArrayRCP<const size_t> host_rowptrs_;
    ArrayRCP<const int>    dev_rowptrs_;
    ArrayRCP<const int>    host_colinds_, dev_colinds_;
    // TODO: add CSC data, for efficient transpose multiply
  };

  /// \class CUSPARSECrsGraph
  /// \brief Sparse matrix class supporting the CUSPARSE library.
  /// \tparam Scalar The type of each entry in the sparse matrix
  /// \tparam Node The Kokkos Node type
  ///
  /// CUSPARSEOps uses this class to represent the local sparse matrix.
  ///
  /// \note This class is really an implementation detail of Tpetra.
  ///   Tpetra users do not normally interact with this class.  In
  ///   fact, Tpetra::CrsMatrix does not give users a way to access an
  ///   instance of this class.
  template <class Scalar, class Node>
  class CUSPARSECrsMatrix : public CrsMatrixBase<Scalar,int,Node> {
  public:
    CUSPARSECrsMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph,
                      const RCP<ParameterList> &params);
    void setValues(const ArrayRCP<const Scalar> &vals);
    void setDeviceData(const ArrayRCP<const Scalar> &devvals);
    inline ArrayRCP<const Scalar> getValues() const;
    inline ArrayRCP<const Scalar> getDevValues() const;
    inline bool isInitialized() const;
    void setAnalyses(const RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                     const RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                     const RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans);
    void getAnalyses(RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
                     RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
                     RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans) const;
  private:
    bool isInitialized_;
    // cusparse analysis handles
    RCP<cusparseSolveAnalysisInfo_t> analysisNoTrans_,
      analysisConjTrans_,
      analysisTrans_;
    // matrix data
    ArrayRCP<const Scalar> vals_, dev_vals_;
    // TODO: add CSC data, for efficient transpose multiply
  };

  template <class Node>
  CUSPARSECrsGraph<Node>::CUSPARSECrsGraph(int numRows, int numCols,
                                           const RCP<Node> &node,
                                           const RCP<ParameterList> &params)
  : CrsGraphBase<int,Node>(numRows,numCols,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    CUSPARSEdetails::Session::init();
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isEmpty() const
  { return isEmpty_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setStructure(const ArrayRCP<const size_t> &ptrs,
                                            const ArrayRCP<const int> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const int numrows = this->getNumRows();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)ptrs.size() != (size_t)numrows+1
        || ptrs[0] != 0
        || (size_t)inds.size() != (size_t)ptrs[numrows],
        std::runtime_error, ": graph data not coherent."
    )
    const int numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, ": matrix has already been initialized"
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    host_rowptrs_ = ptrs;
    host_colinds_ = inds;
    isInitialized_ = true;
  }

  template <class Node>
  ArrayRCP<const size_t> CUSPARSECrsGraph<Node>::getPointers() const
  { return host_rowptrs_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getIndices() const
  { return host_colinds_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getDevPointers() const
  { return dev_rowptrs_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getDevIndices() const
  { return dev_colinds_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setDeviceData(const ArrayRCP<const int> &devptrs,
                                             const ArrayRCP<const int>    &devinds)
  { dev_rowptrs_ = devptrs; dev_colinds_ = devinds; }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isInitialized() const
  { return isInitialized_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    std::string tfecfFuncName("setMatDesc()");
    //
    uplo_ = uplo;
    diag_ = diag;
    //
    matdescr_ = CUSPARSEdetails::createMatDescr();
    cusparseStatus_t stat = cusparseSetMatIndexBase(*matdescr_, CUSPARSE_INDEX_BASE_ZERO);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        stat != CUSPARSE_STATUS_SUCCESS,
        std::runtime_error,
        ": error setting matrix descriptor (index base)."
    )
    // upper or lower
    if (uplo == Teuchos::UPPER_TRI) {
      stat = cusparseSetMatFillMode(*matdescr_, CUSPARSE_FILL_MODE_UPPER);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error,
          ": error setting matrix descriptor (upper)."
      )
    }
    else if (uplo == Teuchos::LOWER_TRI) {
      stat = cusparseSetMatFillMode(*matdescr_, CUSPARSE_FILL_MODE_LOWER);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error,
          ": error setting matrix descriptor (lower)."
      )
    }
    if (diag == Teuchos::UNIT_DIAG) {
      stat = cusparseSetMatDiagType(*matdescr_, CUSPARSE_DIAG_TYPE_UNIT);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error,
          ": error setting matrix descriptor (unit)."
      )
    }
    else {
      stat = cusparseSetMatDiagType(*matdescr_, CUSPARSE_DIAG_TYPE_NON_UNIT);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error,
          ": error setting matrix descriptor (non-unit)."
      )
    }
  }

  template <class Node>
  void CUSPARSECrsGraph<Node>::getMatDesc(Teuchos::EUplo &uplo,
                                          Teuchos::EDiag &diag) const
  {
    uplo = uplo_;
    diag = diag_;
  }

  template <class Node>
  RCP<cusparseMatDescr_t> CUSPARSECrsGraph<Node>::getMatDesc() const
  {
    return matdescr_;
  }

  template <class Scalar, class Node>
  CUSPARSECrsMatrix<Scalar,Node>::CUSPARSECrsMatrix(
                      const RCP<const CUSPARSECrsGraph<Node> > &graph,
                      const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,int,Node>(graph,params)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::setValues(const ArrayRCP<const Scalar> &vals)
  {
    const char tfecfFuncName[] = "setValues";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_ == true, std::runtime_error,
      ": matrix is already initialized.");
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Node>
  ArrayRCP<const Scalar> CUSPARSECrsMatrix<Scalar,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Node>
  bool CUSPARSECrsMatrix<Scalar,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Scalar, class Node>
  ArrayRCP<const Scalar> CUSPARSECrsMatrix<Scalar,Node>::getDevValues() const
  { return dev_vals_; }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::
  setDeviceData(const ArrayRCP<const Scalar> &devvals)
  { dev_vals_ = devvals; }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::
  setAnalyses(const RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
              const RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
              const RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans)
  {
    analysisNoTrans_   = analysisNoTrans;
    analysisTrans_     = analysisTrans;
    analysisConjTrans_ = analysisConjTrans;
  }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::
  getAnalyses(RCP<cusparseSolveAnalysisInfo_t> &analysisNoTrans,
              RCP<cusparseSolveAnalysisInfo_t> &analysisTrans,
              RCP<cusparseSolveAnalysisInfo_t> &analysisConjTrans) const
  {
    analysisNoTrans   = analysisNoTrans_;
    analysisTrans     = analysisTrans_;
    analysisConjTrans = analysisConjTrans_;
  }

  namespace Details {
    /// \class CUSPARSEOpsBase
    /// \brief Base class for cuSPARSE implementation of local sparse operations.
    /// \tparam Node The Kokkos Node type.
    /// \ingroup kokkos_crs_ops
    ///
    /// This class is NOT meant for use as a LocalMatOps
    /// implementation.  It is the base class of CUSPARSEOps and
    /// contains common code for all \c Scalar type specializations of
    /// that class.  Thus, this class is not templated on \c Scalar
    /// and does not know about matrices (which depend on \c Scalar).
    ///
    /// CUSPARSEOps depends on a class hierarchy, because we want
    /// there to be stub implementations for \c Scalar types that the
    /// cuSPARSE library itself does not support.  (This makes it
    /// easier to support explicit template instantiation, by reducing
    /// the number of special cases.)  CUSPARSEOpsBase contains common
    /// code for both stub and actual implementations.  CUSPARSEOps
    /// specializations for unsupported \c Scalar types inherit from
    /// the stub implementation in CUSPARSEOpsStub, and CUSPARSEOps
    /// specializations for \c Scalar types that cuSPARSE supports
    /// inherit from CUSPARSEOpsGeneric.  Since derived classes do not
    /// inherit their base classes' typedefs, we have to put all the
    /// public typedefs in the CUSPARSEOps specializations themselves.
    template <class Node>
    class CUSPARSEOpsBase : public Teuchos::Describable {
    public:
      //@{
      //! @name Typedefs and structs

      //! The type of the (local) indices describing the structure of the sparse matrix.
      typedef int  ordinal_type;
      //! The Kokkos Node type.
      typedef Node node_type;

      /// \brief Typedef for local graph class.
      ///
      /// cuSPARSE does not support arbitrary Ordinal types.  This is
      /// why this class' other_type typedef has a default definition
      /// which results in a compile error if used.  The graph<int,N>
      /// specialization (see below) is defined.
      template <class O, class N>
      struct graph {
        // This will raise a compile error if the given ordinal type O is not supported.
        typedef typename O::this_ordinal_not_supported_by_cusparse graph_type;
      };

      //! Partial specialization of graph<O,N> for O=int.
      template <class N>
      struct graph<int, N> {
        typedef CUSPARSECrsGraph<N> graph_type;
      };

      //@}
      //! @name Constructors/Destructor
      //@{

      //! Constructor accepting and retaining a node object.
      CUSPARSEOpsBase (const RCP<Node> &node) :
        node_ (node),
        numRows_ (0),
        numCols_ (0),
        numNZ_ (0),
        isInitialized_ (false)
      {
        CUSPARSEdetails::Session::init();
        // Check at compile time that Node is a CUDA Node type.
        Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta;
        (void) cta;
      }

      //! Destructor (declared virtual for memory safety of derived classes).
      virtual ~CUSPARSEOpsBase () {}

      //@}
      //! \name Implementation of Teuchos::Describable
      //@{

      //! One-line description of this instance.
      virtual std::string description () const {
        using Teuchos::TypeNameTraits;
        std::ostringstream os;
        os << "KokkosClassic::CUSPARSEOpsBase<"
           << "Node=" << TypeNameTraits<Node>::name()
           << ">";
        return os.str();
      }

      //@}
      //! @name Accessor routines.
      //@{

      //! The Kokkos Node with which this object was instantiated.
      RCP<Node> getNode () const {
        return node_;
      }

      //@}
      //! @name Initialization of graph and matrix
      //@{

      //! Allocate and initialize storage for row offsets.
      static ArrayRCP<size_t>
      allocRowPtrs (const RCP<Node> &node,
                    const ArrayView<const size_t>& numEntriesPerRow)
      {
        // alloc page-locked ("pinned") memory on the host,
        // specially allocated and specially deallocated
        CUDANodeHostPinnedDeallocator<size_t> dealloc;
        ArrayRCP<size_t> ptrs = dealloc.alloc(numEntriesPerRow.size() + 1);
        ptrs[0] = 0;
        std::partial_sum (numEntriesPerRow.getRawPtr(),
                          numEntriesPerRow.getRawPtr() + numEntriesPerRow.size(),
                          ptrs.begin() + 1);
        return ptrs;
      }


      //! Allocate and initialize the storage for a sparse graph.
      template <class T>
      static ArrayRCP<T>
      allocStorage (const RCP<Node> &node,
                    const ArrayView<const size_t>& rowPtrs)
      {
        // alloc page-locked ("pinned") memory on the host,
        // specially allocated and specially deallocated
        const int totalNumEntries = *(rowPtrs.end()-1);
        CUDANodeHostPinnedDeallocator<T> dealloc;
        ArrayRCP<T> buf = dealloc.alloc (totalNumEntries);
        std::fill (buf.begin(), buf.end(), Teuchos::ScalarTraits<T>::zero());
        return buf;
      }

      //! Finalize a graph is null for CUSPARSE.
      static void
      finalizeGraph (Teuchos::EUplo uplo, Teuchos::EDiag diag,
                     CUSPARSECrsGraph<Node> &graph,
                     const RCP<ParameterList> &params)
      {
        const size_t CUDA_MAX_INT = 2147483647;
        const char prefix[] = "finalizeGraph";
        RCP<Node> node = graph.getNode ();
        TEUCHOS_TEST_FOR_EXCEPTION(
           graph.isInitialized() == false, std::runtime_error,
           prefix << ": graph has not yet been initialized.");

        // diag: have to allocate and indicate
        ArrayRCP<int> devinds, devptrs;
        ArrayRCP<const int> hostinds = graph.getIndices ();
        ArrayRCP<const size_t> hostptrs = graph.getPointers ();
        const int numRows = graph.getNumRows ();
        // set description
        graph.setMatDesc (uplo, diag);
        // allocate and initialize data
        if (diag == Teuchos::UNIT_DIAG) {
          // CUSPARSE, unfortunately, always assumes that the diagonal
          // entries are present in the storage.  Therefore, this flag
          // only specifies whether they are considered or not; they
          // are assumed to be present, and neglecting them will
          // result in incorrect behavior (causing the adjacent entry
          // to be neglected instead).  However, our API doesn't give
          // us explicit diagonal entries if diag ==
          // Teuchos::UNIT_DIAG; this adaptor must therefore allocate
          // space for them.  Furthermore, because there is no support
          // in our API or CUSPARSE to ignore them on multiply, we
          // must set the values of the explicit diagonals to zero
          // (making them effectively ignored on multiply)
          const size_t numnz = hostinds.size() + numRows;
          TEUCHOS_TEST_FOR_EXCEPTION(
            numnz > CUDA_MAX_INT, std::runtime_error,
            "KokkosClassic::CUSPARSEOps: CUSPARSE does not support more than "
            << CUDA_MAX_INT << " non-zeros.");

          devptrs = node->template allocBuffer<int>( numRows+1 );
          if (numnz) {
            devinds = node->template allocBuffer<int>( numnz );
          }
          ArrayRCP<int> h_devptrs =
            node->viewBufferNonConst (WriteOnly, numRows+1, devptrs);
          ArrayRCP<int> h_devinds =
            node->viewBufferNonConst (WriteOnly, numnz, devinds);
          for (int r = 0; r < numRows; ++r) {
            h_devptrs[r] = static_cast<int> (hostptrs[r] + r);
            if (uplo == Teuchos::LOWER_TRI) {
              // copy the explicit entries, then set the last one manually
              std::copy (hostinds.begin()+hostptrs[r],
                         hostinds.begin()+hostptrs[r+1],
                         h_devinds.begin()+h_devptrs[r]);
              h_devinds[hostptrs[r+1]+r+1-1] = r;
            }
            else {
              // set the first entry, then skip it in the copy
              h_devinds[h_devptrs[r]] = r;
              std::copy (hostinds.begin()+hostptrs[r],
                         hostinds.begin()+hostptrs[r+1],
                         h_devinds.begin()+h_devptrs[r]+1);
            }
          }
          h_devptrs[numRows] = static_cast<int> (hostptrs[numRows] + numRows);
          // copy back
          h_devptrs = null;
          h_devinds = null;
        }
        else {
          // our format == their format; just allocate and copy
          const size_t numnz = hostinds.size();
          TEUCHOS_TEST_FOR_EXCEPTION(
            numnz > CUDA_MAX_INT, std::runtime_error,
            "KokkosClassic::CUSPARSEOps: CUSPARSE does not support more than "
            << CUDA_MAX_INT << " non-zeros.");

          devptrs = node->template allocBuffer<int> (numRows + 1);
          ArrayRCP<int> h_devptrs =
            node->viewBufferNonConst (WriteOnly, numRows+1, devptrs);
          std::copy (hostptrs.begin(), hostptrs.end(), h_devptrs.begin());
          h_devptrs = null;
          if (numnz) {
            devinds = node->template allocBuffer<int> (numnz);
            node->copyToBuffer (numnz, hostinds (), devinds);
          }
        }
        // set the data
        graph.setDeviceData (devptrs, devinds);
      }
      //@}

    private:
      //! Copy constructor (protected and unimplemented)
      CUSPARSEOpsBase (const CUSPARSEOpsBase& source);

      //! The Kokkos Node instance given to this object's constructor.
      RCP<Node> node_;

    protected:
      ArrayRCP<const int> rowPtrs_;
      ArrayRCP<const int> colInds_;
      int numRows_;
      int numCols_;
      int numNZ_;
      bool isInitialized_;
    };

    /// \class CUSPARSEOpsStub
    /// \brief Stub version of cuSPARSE implementation of local sparse operations.
    /// \tparam Scalar The type of entries of the sparse matrix.
    /// \tparam Node The Kokkos Node type.
    /// \ingroup kokkos_crs_ops
    ///
    /// This class is helpful for avoiding link errors when doing
    /// explicit instantiation for Scalar types that cuSPARSE does not
    /// support.  Currently, cuSPARSE only supports <tt>float</tt> and
    /// <tt>double</tt>.
    template <class Scalar, class Node>
    class CUSPARSEOpsStub : public Details::CUSPARSEOpsBase<Node> {
    public:

      typedef int Ordinal;

      //! @name Constructors/Destructor
      //@{

      //! Constructor accepting and retaining a node object.
      CUSPARSEOpsStub (const RCP<Node> &node) :
        CUSPARSEOpsBase<Node> (node) {}

      /// \brief "Sum constructor": compute *this = alpha*A + beta*B.
      ///
      /// The resulting matrix shares the Node instance and copies the
      /// parameters of the matrix A.
      CUSPARSEOpsStub (const Scalar& alpha,
                       const CUSPARSEOpsStub<Scalar, Node>& A,
                       const Scalar& beta,
                       const CUSPARSEOpsStub<Scalar, Node>& B)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "CUSPARSEOps: sum constructor not implemented.");
      }

      //! Destructor
      virtual ~CUSPARSEOpsStub () {}

      //@}
      //! \name Implementation of Teuchos::Describable
      //@{

      //! One-line description of this instance.
      virtual std::string description () const {
        using Teuchos::TypeNameTraits;
        std::ostringstream os;
        os << "KokkosClassic::CUSPARSEOps<"
           << "Scalar=" << TypeNameTraits<Scalar>::name()
           << ", Node=" << TypeNameTraits<Node>::name()
           << ">";
        return os.str();
      }

      //@}
      /// \name Initialization of graph and matrix.
      ///
      /// Many methods in this category are inherited from the base class.
      //@{

      static void
      finalizeGraph (Teuchos::EUplo uplo, Teuchos::EDiag diag,
                     CUSPARSECrsGraph<Node> &graph,
                     const RCP<ParameterList> &params)
      {
        CUSPARSEOpsBase<Node>::finalizeGraph (uplo, diag, graph, params);
      }

      //! Finalize the matrix of an already-finalized graph.
      static void
      finalizeMatrix (const CUSPARSECrsGraph<Node> &graph,
                      CUSPARSECrsMatrix<Scalar,Node> &matrix,
                      const RCP<ParameterList> &params)
      {
        (void) graph;
        (void) matrix;
        (void) params;
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "KokkosClassic::CUSPARSEOps::"
          "finalizeMatrix: Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      //! Finalize a graph and a matrix.
      static void
      finalizeGraphAndMatrix (Teuchos::EUplo uplo, Teuchos::EDiag diag,
                              CUSPARSECrsGraph<Node> &graph,
                              CUSPARSECrsMatrix<Scalar,Node> &matrix,
                              const RCP<ParameterList> &params)
      {
        (void) diag;
        (void) graph;
        (void) matrix;
        (void) params;
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "KokkosClassic::CUSPARSEOps::"
          "finalizeGraphAndMatrix: Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      //! Initialize sparse operations with a graph and matrix
      void
      setGraphAndMatrix (const RCP<const CUSPARSECrsGraph<Node> > &graph,
                         const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix)
      {
        (void) graph;
        (void) matrix;
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "KokkosClassic::CUSPARSEOps::"
          "setGraphAndMatrix: Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      //@}
      //! @name Computational methods
      //@{

      /// \brief Y := alpha * Op(A) * X.
      ///
      /// Apply the local sparse matrix A (or its transpose or conjugate
      /// transpose) to a multivector X, overwriting Y with the result.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to apply the matrix, its transpose,
      ///   or its conjugate transpose (if applicable).
      ///
      /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
      ///   multiply the result: \f$Y := \alpha A X\f$.
      ///
      /// \param X [in] Input multivector.
      ///
      /// \param Y [out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      multiply (Teuchos::ETransp trans,
                RangeScalar alpha,
                const MultiVector<DomainScalar,Node> &X,
                MultiVector<RangeScalar,Node> &Y) const
      {
        (void) trans;
        (void) alpha;
        (void) X;
        (void) Y;
        const char tfecfFuncName[] = "KokkosClassic::CUSPARSEOps::multiply";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::runtime_error, ": Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      /// \brief Y := Y + alpha * Op(A) * X.
      ///
      /// Apply the local sparse matrix A (or its transpose or conjugate
      /// transpose) to a multivector X, accumulating the result into Y.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to apply the matrix, its transpose,
      ///   or its conjugate transpose (if applicable).
      ///
      /// \param alpha [in] Scalar coefficient \f$\alpha\f$ on the
      ///   product.
      ///
      /// \param X [in] Input multivector.
      ///
      /// \param beta [in] Scalar coefficient \f$\beta\f$ on the
      ///   accumulation in \f$Y\f$.
      ///
      /// \param Y [in/out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      multiply (Teuchos::ETransp trans,
                RangeScalar alpha,
                const MultiVector<DomainScalar,Node> &X,
                RangeScalar beta,
                MultiVector<RangeScalar,Node> &Y) const
      {
        (void) trans;
        (void) alpha;
        (void) X;
        (void) beta;
        (void) Y;
        const char tfecfFuncName[] = "KokkosClassic::CUSPARSEOps::multiply";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::runtime_error, ": Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
      ///
      /// Solve the (upper or lower) triangular system Y = Op(A) X.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to solve with the matrix, its
      ///   transpose, or its conjugate transpose (if applicable).
      ///
      /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
      ///   else LOWER_TRI if the matrix is lower triangular.
      ///
      /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
      ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
      ///
      /// \param Y [in] Input multivector.
      ///
      /// \param X [out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      solve (Teuchos::ETransp trans,
             const MultiVector<DomainScalar,Node> &Y,
             MultiVector<RangeScalar,Node> &X) const
      {
        (void) trans;
        (void) Y;
        (void) X;
        const char tfecfFuncName[] = "KokkosClassic::CUSPARSEOps::solve";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::runtime_error, ": Not implemented for Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << ".");
      }

      template <class DomainScalar, class RangeScalar>
      void
      gaussSeidel (const MultiVector<DomainScalar,Node> &B,
                   MultiVector< RangeScalar,Node> &X,
                   const MultiVector<Scalar,Node> &D,
                   const RangeScalar& dampingFactor,
                   const ESweepDirection direction) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "CUSPARSEOps: gaussSeidel not implemented");
      }

      template <class DomainScalar, class RangeScalar>
      void
      reorderedGaussSeidel (const MultiVector<DomainScalar,Node> &B,
			    MultiVector< RangeScalar,Node> &X,
			    const MultiVector<Scalar,Node> &D,
			    const ArrayView<Ordinal> & rowIndices,
			    const RangeScalar& dampingFactor,
			    const ESweepDirection direction) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "CUSPARSEOps: reorderedGaussSeidel not implemented");
      }


      /// \brief "Add in place": compute <tt>*this = alpha*A + beta*(*this)</tt>.
      ///
      /// This method may choose to reuse storage of <tt>*this</tt>.
      void
      addInPlace (const Scalar& alpha,
                  const CUSPARSEOpsStub<Scalar, Node>& A,
                  const Scalar& beta)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "CUSPARSEOps: addInPlace not implemented");
      }
      //@}

    private:
      //! Copy constructor (protected and unimplemented)
      CUSPARSEOpsStub (const CUSPARSEOpsStub& source);
    };

    /// \class CUSPARSEOpsGeneric
    /// \brief "Generic" cuSPARSE implementation of local sparse operations.
    /// \tparam Scalar The type of entries of the sparse matrix.
    /// \tparam Node The Kokkos Node type.
    /// \ingroup kokkos_crs_ops
    ///
    /// This is the "generic" cuSPARSE implementation of local sparse
    /// operations.  We put "generic" in quotes because cuSPARSE
    /// itself only supports \c float and \c double data types.  This
    /// class wraps both of these implementations.  The actual
    /// CUSPARSEOps class inherits from this class.
    ///
    /// This class is one of various classes in Kokkos that implement
    /// local sparse matrix-(multi)vector multiply and sparse
    /// triangular solve.  ("Local" means "on a single node; not using
    /// MPI.")  Examples include DefaultHostSparseOps, AltSparseOps,
    /// and MklSparseOps for host-based Kokkos Nodes, and this class
    /// and CuspOps for NVIDIA GPUs (Graphics Processing Units).  This
    /// class provides an interface to the local sparse operations
    /// provided by <a
    /// href="http://developer.nvidia.com/cuda/cusparse">cuSPARSE</a>,
    /// the NVIDIA CUDA Sparse Matrix library.
    ///
    /// \note Unlike CuspOps and the other local sparse operations
    ///   classes mentioned above, this class is not templated on the
    ///   Ordinal type (of column indices in the sparse matrix).  This
    ///   is because cuSPARSE currently only supports column indices of
    ///   type \c int.
    template <class Scalar, class Node>
    class CUSPARSEOpsGeneric : public Details::CUSPARSEOpsBase<Node> {
    public:

      typedef int Ordinal;

      //! @name Constructors/Destructor
      //@{

      //! Constructor accepting and retaining a node object.
      CUSPARSEOpsGeneric (const RCP<Node> &node) :
        CUSPARSEOpsBase<Node> (node) {}

      /// \brief "Sum constructor": compute *this = alpha*A + beta*B.
      ///
      /// The resulting matrix shares the Node instance and copies the
      /// parameters of the matrix A.
      CUSPARSEOpsGeneric (const Scalar& alpha,
                          const CUSPARSEOpsGeneric<Scalar, Node>& A,
                          const Scalar& beta,
                          const CUSPARSEOpsGeneric<Scalar, Node>& B)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "CUSPARSEOps: sum constructor not implemented.");
      }

      //! Destructor
      virtual ~CUSPARSEOpsGeneric () {}

      //@}
      //! \name Implementation of Teuchos::Describable
      //@{

      //! One-line description of this instance.
      virtual std::string description () const {
        using Teuchos::TypeNameTraits;
        std::ostringstream os;
        os << "KokkosClassic::CUSPARSEOpsGeneric<"
           << "Scalar=" << TypeNameTraits<Scalar>::name()
           << ", Node=" << TypeNameTraits<Node>::name()
           << ">";
        return os.str();
      }

      //@}
      /// \name Initialization of graph and matrix.
      ///
      /// Many methods in this category are inherited from the base class.
      //@{

      static void
      finalizeGraph (Teuchos::EUplo uplo, Teuchos::EDiag diag,
                     CUSPARSECrsGraph<Node> &graph,
                     const RCP<ParameterList> &params)
      {
        CUSPARSEOpsBase<Node>::finalizeGraph (uplo, diag, graph, params);
      }

      //! Finalize the matrix of an already-finalized graph.
      static void
      finalizeMatrix (const CUSPARSECrsGraph<Node> &graph,
                      CUSPARSECrsMatrix<Scalar,Node> &matrix,
                      const RCP<ParameterList> &params)
      {
        RCP<Node> node = graph.getNode();
        TEUCHOS_TEST_FOR_EXCEPTION(
           matrix.isInitialized() == false, std::runtime_error,
           "KokkosClassic::CUSPARSEOps::finalizeMatrix: "
           "matrix has not yet been initialized.");

        // diag: have to allocate and indicate
        Teuchos::EUplo uplo;
        Teuchos::EDiag diag;
        graph.getMatDesc(uplo,diag);
        ArrayRCP<Scalar>       devvals;
        ArrayRCP<const size_t> hostptrs = graph.getPointers();
        ArrayRCP<const Scalar> hostvals = matrix.getValues();
        const int numRows = graph.getNumRows();
        if (diag == Teuchos::UNIT_DIAG) {
          // CUSPARSE requires explicit diagonals, even if UNIT_DIAG
          // our API requires no explicit diagonals when UNIT_DIAG
          // therefore, we have to add them during our copy to device.
          // they must be set to zero, however, to prevent them from
          // affecting the result of a multiplication (which ignore UNIT_DIAG)
          const int numnz = hostptrs[numRows] + numRows;
          devvals = node->template allocBuffer<Scalar>( numnz );
          ArrayRCP<Scalar> h_devvals =
            node->viewBufferNonConst (WriteOnly, numnz, devvals);
          for (int r = 0; r < numRows; ++r) {
            if (uplo == Teuchos::LOWER_TRI) {
              std::copy (hostvals.begin()+hostptrs[r], hostvals.begin()+hostptrs[r+1],
                         h_devvals.begin()+hostptrs[r]+r );
              // diag goes at end
              h_devvals[hostptrs[r+1]+r+1-1] = Teuchos::ScalarTraits<Scalar>::zero();
            }
            else {
              // diag goes at beginning
              h_devvals[hostptrs[r]+r] = Teuchos::ScalarTraits<Scalar>::zero();
              std::copy( hostvals.begin()+hostptrs[r], hostvals.begin()+hostptrs[r+1],
                         h_devvals.begin()+hostptrs[r]+r+1 );
            }
          }
          // copy back
          h_devvals = null;
        }
        else {
          const int numnz = hostptrs[numRows];
          if (numnz) {
            devvals = node->template allocBuffer<Scalar>( numnz );
            node->copyToBuffer (numnz, hostvals (), devvals);
          }
        }
        matrix.setDeviceData (devvals);

        const int numnz = static_cast<int> (devvals.size ());
        RCP<cusparseMatDescr_t> descr = graph.getMatDesc ();
        ArrayRCP<const int> devptrs = graph.getDevPointers (),
          devinds = graph.getDevIndices ();

        RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
        // look at the parameter list and do any analyses requested for solves
        RCP<cusparseSolveAnalysisInfo_t> ai_non, ai_trans, ai_conj;
        if (params != null && params->get("Prepare Solve",false)) {
          ai_non = CUSPARSEdetails::createSolveAnalysisInfo();
          cusparseStatus_t stat =
            CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(*hndl, CUSPARSE_OPERATION_NON_TRANSPOSE,numRows, numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(), *ai_non);

          TEUCHOS_TEST_FOR_EXCEPTION(
            stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
            "KokkosClassic::CUSPARSEOps::finalizeMatrix: "
            "CSRSM_analysis(non-trans) returned error " << stat);
        }
        if (params != null && params->get ("Prepare Transpose Solve", false)) {
          ai_trans = CUSPARSEdetails::createSolveAnalysisInfo();
          cusparseStatus_t stat =
            CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(
                                                                              *hndl, CUSPARSE_OPERATION_TRANSPOSE,numRows,
                                                                              numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(),
                                                                              *ai_trans
                                                                              );
          TEUCHOS_TEST_FOR_EXCEPTION(
            stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
            "KokkosClassic::CUSPARSEOps::finalizeMatrix: "
            "CSRSM_analysis(trans) returned error " << stat);
        }
        if (params != null && params->get ("Prepare Conjugate Transpose Solve", false)) {
          ai_conj = CUSPARSEdetails::createSolveAnalysisInfo ();
          cusparseStatus_t stat =
            CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_analysis(
                                                                              *hndl, CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE,numRows,
                                                                              numnz, *descr, devvals.getRawPtr(), devptrs.getRawPtr(), devinds.getRawPtr(),
                                                                              *ai_conj
                                                                              );
          TEUCHOS_TEST_FOR_EXCEPTION(
            stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
            "KokkosClassic::CUSPARSEOps::finalizeMatrix: "
            "CSRSM_analysis(conj-trans) returned error " << stat);
        }
        matrix.setAnalyses (ai_non, ai_trans, ai_conj);
        //
        // if (params != null && params->get("Prepare Transpose Multiply",false)) {
        //   // finish: compute CSC for transpose
        //   TEUCHOS_TEST_FOR_EXCEPT(true);
        // }
      }

      //! Finalize a graph and a matrix.
      static void
      finalizeGraphAndMatrix (Teuchos::EUplo uplo, Teuchos::EDiag diag,
                              CUSPARSECrsGraph<Node> &graph,
                              CUSPARSECrsMatrix<Scalar,Node> &matrix,
                              const RCP<ParameterList> &params)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          graph.isInitialized() == false, std::runtime_error,
          "KokkosClassic::CUSPARSEOps::finalizeGraphAndMatrix: "
          "graph has not yet been initialized.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          matrix.isInitialized() == false, std::runtime_error,
          "KokkosClassic::CUSPARSEOps::finalizeGraphAndMatrix: "
          "matrix has not yet been initialized.");
        // no benefit to doing them together; do them separately
        finalizeGraph (uplo, diag, graph, params);
        finalizeMatrix (graph, matrix, params);
      }

      //! Initialize sparse operations with a graph and matrix
      void
      setGraphAndMatrix (const RCP<const CUSPARSECrsGraph<Node> > &graph,
                         const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix)
      {
        const char tfecfFuncName[] = "setGraphAndMatrix";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          this->isInitialized_ == true, std::runtime_error,
          ": operators are already initialized.");
        // Get cuSPARSE data from the matrix
        this->numRows_ = graph->getNumRows ();
        this->numCols_ = graph->getNumCols ();
        this->matdescr_ = graph->getMatDesc ();
        this->rowPtrs_ = graph->getDevPointers ();
        this->colInds_ = graph->getDevIndices ();
        this->rowVals_ = matrix->getDevValues ();
        this->numNZ_ = this->colInds_.size ();
        matrix->getAnalyses (aiNoTrans_, aiTrans_, aiConjTrans_);
        this->isInitialized_ = true;
      }

      //@}
      //! @name Computational methods
      //@{

      /// \brief Y := alpha * Op(A) * X.
      ///
      /// Apply the local sparse matrix A (or its transpose or conjugate
      /// transpose) to a multivector X, overwriting Y with the result.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to apply the matrix, its transpose,
      ///   or its conjugate transpose (if applicable).
      ///
      /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
      ///   multiply the result: \f$Y := \alpha A X\f$.
      ///
      /// \param X [in] Input multivector.
      ///
      /// \param Y [out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      multiply (Teuchos::ETransp trans,
                RangeScalar alpha,
                const MultiVector<DomainScalar,Node> &X,
                MultiVector<RangeScalar,Node> &Y) const
      {
        // CUSPARSE doesn't support mixed precision
        Teuchos::CompileTimeAssert<
          Teuchos::TypeTraits::is_same<DomainScalar,Scalar>::value == false ||
          Teuchos::TypeTraits::is_same< RangeScalar,Scalar>::value == false> cta;
        (void) cta;

        const char tfecfFuncName[] = "multiply(trans,alpha,X,Y)";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          this->isInitialized_ == false, std::runtime_error,
          ": sparse operators have not been initialized with graph and matrix data;"
          "call setGraphAndMatrix() first.");
        // get pointers,stride from X and Y
        int stride_x = (int)X.getStride(),
          stride_y = (int)Y.getStride();
        const Scalar * data_x = X.getValues().getRawPtr();
        Scalar * data_y = Y.getValuesNonConst().getRawPtr();
        const int numMatRows = this->numRows_;
        const int numMatCols = this->numCols_;
        const int opRows     = (trans == Teuchos::NO_TRANS ? numMatRows : numMatCols);
        const int opCols     = (trans == Teuchos::NO_TRANS ? numMatCols : numMatRows);
        const int numRHS     = X.getNumCols();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          X.getNumCols() != Y.getNumCols(), std::runtime_error,
          ": X and Y do not have the same number of column vectors.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)X.getNumRows() != (size_t)opCols, std::runtime_error,
          ": size of X is not congruous with dimensions of operator.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)Y.getNumRows() != (size_t)opRows, std::runtime_error,
          ": size of Y is not congruous with dimensions of operator.");
        // CUSPARSE will short-circuit on these without doing anything,
        // so we have to do it
        if (opRows == 0.0 || opCols == 0.0) {
          // Y <- alpha*A*X == alpha*0*X == 0
          KokkosClassic::DefaultArithmetic< MultiVector<DomainScalar,Node> >
            ::Init(Y, Teuchos::ScalarTraits<Scalar>::zero() );
          return;
        }
        // call mat-vec
        cusparseOperation_t op;
        if      (trans == Teuchos::NO_TRANS)     op = CUSPARSE_OPERATION_NON_TRANSPOSE;
        else if (trans == Teuchos::TRANS)        op = CUSPARSE_OPERATION_TRANSPOSE;
        else /* (trans == Teuchos::CONJ_TRANS)*/ op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        RCP<const cusparseHandle_t> sess = CUSPARSEdetails::Session::getHandle();
        const Scalar s_alpha = (Scalar)alpha,
          s_beta  = Teuchos::ScalarTraits<Scalar>::zero();
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        cudaError_t err = cudaGetLastError();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaGetLastError() returned error after function call:\n"
          << cudaGetErrorString (err));
        err = cudaThreadSynchronize();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaThreadSynchronize() returned error after function call:\n"
          << cudaGetErrorString(err));
#endif // HAVE_KOKKOSCLASSIC_DEBUG
        cusparseStatus_t stat;
        // we're only ever multiplying by a general matrix
        stat = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_GENERAL);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error, ": error setting matrix descriptor (general).");
        // do the multiply
        stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRMM(
          *sess, op, numMatRows, numRHS, numMatCols, this->numNZ_, &s_alpha,
          *matdescr_, this->rowVals_.getRawPtr(), this->rowPtrs_.getRawPtr(), this->colInds_.getRawPtr(),
          data_x, stride_x, &s_beta, data_y, stride_y);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        err = cudaGetLastError();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaGetLastError() returned error after function call:\n"
          << cudaGetErrorString(err) );
        err = cudaThreadSynchronize();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaThreadSynchronize() returned error after function call:\n"
          << cudaGetErrorString(err) );
#endif
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error, ": CSRMM returned error " << stat);
        return;
      }

      /// \brief Y := Y + alpha * Op(A) * X.
      ///
      /// Apply the local sparse matrix A (or its transpose or conjugate
      /// transpose) to a multivector X, accumulating the result into Y.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to apply the matrix, its transpose,
      ///   or its conjugate transpose (if applicable).
      ///
      /// \param alpha [in] Scalar coefficient \f$\alpha\f$ on the
      ///   product.
      ///
      /// \param X [in] Input multivector.
      ///
      /// \param beta [in] Scalar coefficient \f$\beta\f$ on the
      ///   accumulation in \f$Y\f$.
      ///
      /// \param Y [in/out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      multiply (Teuchos::ETransp trans,
                RangeScalar alpha,
                const MultiVector<DomainScalar,Node> &X,
                RangeScalar beta,
                MultiVector<RangeScalar,Node> &Y) const
      {
        using Teuchos::TypeTraits::is_same;

        // CUSPARSE doesn't support mixed precision
        Teuchos::CompileTimeAssert<is_same<DomainScalar,Scalar>::value == false ||
                                   is_same<RangeScalar,Scalar>::value == false> cta;
        (void) cta;

        const char tfecfFuncName[] = "multiply(trans,alpha,X,beta,Y)";
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          this->isInitialized_ == false, std::runtime_error,
          ": sparse operators have not been initialized with graph and matrix data;"
          " call setGraphAndMatrix() first.");

        // get pointers,stride from X and Y
        int stride_x = (int)X.getStride(),
          stride_y = (int)Y.getStride();
        const Scalar * data_x = X.getValues().getRawPtr();
        Scalar * data_y = Y.getValuesNonConst().getRawPtr();
        const int numMatRows = this->numRows_;
        const int numMatCols = this->numCols_;
        const int opRows     = (trans == Teuchos::NO_TRANS ? numMatRows : numMatCols);
        const int opCols     = (trans == Teuchos::NO_TRANS ? numMatCols : numMatRows);
        const int numRHS     = X.getNumCols();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          X.getNumCols() != Y.getNumCols(),
          std::runtime_error, ": X and Y do not have the same number of column vectors.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)X.getNumRows() != (size_t)opCols,
          std::runtime_error, ": size of X is not congruous with dimensions of operator.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)Y.getNumRows() != (size_t)opRows,
          std::runtime_error, ": size of Y is not congruous with dimensions of operator.");
        // CUSPARSE will short-circuit on these without doing anything, so we have to do it
        if (opRows == 0.0 || opCols == 0.0) {
          // Y <- alpha*A*X + beta*Y == alpha*0*X + beta*Y == beta*Y
          KokkosClassic::DefaultArithmetic<MultiVector<DomainScalar, Node> >::Scale (Y, beta);
          return;
        }
        // call mat-vec
        cusparseOperation_t op;
        RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
        if      (trans == Teuchos::NO_TRANS)     op = CUSPARSE_OPERATION_NON_TRANSPOSE;
        else if (trans == Teuchos::TRANS)        op = CUSPARSE_OPERATION_TRANSPOSE;
        else /* (trans == Teuchos::CONJ_TRANS)*/ op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        const Scalar s_alpha = (Scalar)alpha,
          s_beta  = (Scalar)beta;
        cusparseStatus_t stat;
        // we're only ever multiplying by a general matrix
        stat = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_GENERAL);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error, ": error setting matrix descriptor (general).");

        stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRMM(
          *hndl, op, numMatRows, numRHS, numMatCols, this->numNZ_, &s_alpha,
          *matdescr_, this->rowVals_.getRawPtr(), this->rowPtrs_.getRawPtr(), this->colInds_.getRawPtr(),
          data_x, stride_x, &s_beta, data_y, stride_y);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        cudaError_t err = cudaGetLastError();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaGetLastError() returned error after function call:\n"
          << cudaGetErrorString(err) );
        err = cudaThreadSynchronize();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          cudaSuccess != err, std::runtime_error,
          ": cudaThreadSynchronize() returned error after function call:\n"
          << cudaGetErrorString(err) );
#endif
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error, ": CSRMM returned error " << stat);
        return;
      }

      /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
      ///
      /// Solve the (upper or lower) triangular system Y = Op(A) X.
      /// Op(A) means A, the transpose of A, or the conjugate transpose
      /// of A, depending on the \c trans argument.
      ///
      /// \tparam DomainScalar The type of entries in the input
      ///   multivector X.  This may differ from the type of entries in
      ///   A or in Y.
      ///
      /// \tparam RangeScalar The type of entries in the output
      ///   multivector Y.  This may differ from the type of entries in
      ///   A or in X.
      ///
      /// \param trans [in] Whether to solve with the matrix, its
      ///   transpose, or its conjugate transpose (if applicable).
      ///
      /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
      ///   else LOWER_TRI if the matrix is lower triangular.
      ///
      /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
      ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
      ///
      /// \param Y [in] Input multivector.
      ///
      /// \param X [out] Result multivector.
      template <class DomainScalar, class RangeScalar>
      void
      solve (Teuchos::ETransp trans,
             const MultiVector<DomainScalar,Node> &Y,
             MultiVector<RangeScalar,Node> &X) const
      {
        using Teuchos::TypeTraits::is_same;
        const char tfecfFuncName[] = "solve";

        // CUSPARSE doesn't support mixed precision;
        // partial specialize, then nix the generic versions
        Teuchos::CompileTimeAssert<is_same<DomainScalar,Scalar>::value == false ||
                                   is_same< RangeScalar,Scalar>::value == false> cta;
        (void) cta;

        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          this->isInitialized_ == false, std::runtime_error,
          ": sparse operators have not been initialized with graph and matrix data;"
          " call setGraphAndMatrix() first.");

        // get pointers,stride from X and Y
        int stride_x = (int)X.getStride(),
          stride_y = (int)Y.getStride();
        const Scalar * data_y = Y.getValues().getRawPtr();
        Scalar * data_x = X.getValuesNonConst().getRawPtr();
        const int numMatRows = X.getNumRows();
        const int numRHS     = X.getNumCols();
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          X.getNumRows() != Y.getNumRows(),
          std::runtime_error, ": X and Y do not have the same number of row vectors.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          X.getNumCols() != Y.getNumCols(),
          std::runtime_error, ": X and Y do not have the same number of column vectors.");
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          Y.getNumRows() != static_cast<size_t> (this->numRows_),
          std::runtime_error,
          ": Y does not have the same number of rows as does the matrix.");

        RCP<cusparseSolveAnalysisInfo_t> solveInfo;
        cusparseOperation_t                     op;
        if (trans == Teuchos::NO_TRANS) {
          solveInfo = aiNoTrans_;
          op = CUSPARSE_OPERATION_NON_TRANSPOSE;
        }
        else if (trans == Teuchos::TRANS) {
          solveInfo = aiTrans_;
          op = CUSPARSE_OPERATION_TRANSPOSE;
        }
        else if (trans == Teuchos::CONJ_TRANS) {
          solveInfo = aiConjTrans_;
          op = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            true, std::runtime_error,
            ": invalid transformation (must be one of NO_TRANS, TRANS or CONJ_TRANS): "
            << trans);
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          solveInfo == null, std::runtime_error,
          ": solve info not computed at matrix finalization time for requested "
          "transformation: " << trans);
        RCP<const cusparseHandle_t> hndl = CUSPARSEdetails::Session::getHandle();
        const Scalar s_alpha = Teuchos::ScalarTraits<Scalar>::one();
        cusparseStatus_t stat;
        // we're only ever solving against a triangular matrix
        stat = cusparseSetMatType(*matdescr_, CUSPARSE_MATRIX_TYPE_TRIANGULAR);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
          ": error setting matrix descriptor (triangular).");

        stat = CUSPARSEdetails::CUSPARSETemplateAdaptors<Scalar>::CSRSM_solve(
          *hndl, op, numMatRows, numRHS, &s_alpha, *matdescr_,
          this->rowVals_.getRawPtr(), this->rowPtrs_.getRawPtr(), this->colInds_.getRawPtr(),
          *solveInfo, data_y, stride_y, data_x, stride_x);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          stat != CUSPARSE_STATUS_SUCCESS,
          std::runtime_error, ": CSRSM_solve returned error " << stat);
        return;
      }

      template <class DomainScalar, class RangeScalar>
      void
      gaussSeidel (const MultiVector<DomainScalar,Node> &B,
                   MultiVector< RangeScalar,Node> &X,
                   const MultiVector<Scalar,Node> &D,
                   const RangeScalar& dampingFactor,
                   const ESweepDirection direction) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "KokkosClassic::CUSPARSEOps::gaussSeidel: Not implemented");
      }

      template <class DomainScalar, class RangeScalar>
      void
      reorderedGaussSeidel (const MultiVector<DomainScalar,Node> &B,
			    MultiVector< RangeScalar,Node> &X,
			    const MultiVector<Scalar,Node> &D,
			    const ArrayView<Ordinal> & rowIndices,
			    const RangeScalar& dampingFactor,
			    const ESweepDirection direction) const
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "KokkosClassic::CUSPARSEOps::reorderedGaussSeidel:  Not implemented");
      }

      /// \brief "Add in place": compute <tt>*this = alpha*A + beta*(*this)</tt>.
      ///
      /// This method may choose to reuse storage of <tt>*this</tt>.
      void
      addInPlace (const Scalar& alpha,
                  const CUSPARSEOpsGeneric<Scalar, Node>& A,
                  const Scalar& beta)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "KokkosClassic::CUSPARSEOps::addinPlace: Not implemented");
      }
      //@}

    private:
      //! Copy constructor (protected and unimplemented)
      CUSPARSEOpsGeneric (const CUSPARSEOpsGeneric& source);

      ArrayRCP<const Scalar> rowVals_;
      RCP<cusparseMatDescr_t> matdescr_;
      RCP<cusparseSolveAnalysisInfo_t> aiNoTrans_, aiTrans_, aiConjTrans_;
    };
  } // namespace Details


  /// \class CUSPARSEOps
  /// \brief cuSPARSE implementation of local sparse operations.
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  /// \ingroup kokkos_crs_ops
  ///
  /// This class is one of various classes in Kokkos that implement
  /// local sparse matrix-(multi)vector multiply and sparse triangular
  /// solve.  ("Local" means "on a single node; not using MPI.")
  /// Examples include DefaultHostSparseOps, AltSparseOps, and
  /// MklSparseOps for host-based Kokkos Nodes, and this class and
  /// CuspOps for NVIDIA GPUs (Graphics Processing Units).  This class
  /// provides an interface to the local sparse operations provided by
  /// <a href="http://developer.nvidia.com/cuda/cusparse">cuSPARSE</a>,
  /// the NVIDIA CUDA Sparse Matrix library.
  ///
  /// \note Unlike CuspOps and the other local sparse operations
  ///   classes mentioned above, this class is not templated on the
  ///   Ordinal type (of column indices in the sparse matrix).  This
  ///   is because cuSPARSE currently only supports column indices of
  ///   type \c int.
  template <class Scalar, class Node>
  class CUSPARSEOps : public Details::CUSPARSEOpsStub<Scalar, Node> {
  public:
    //@{
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef int     ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef CUSPARSEOps<Scalar,Node> sparse_ops_type;

    /// \brief Typedef for local graph class.
    ///
    /// CUSPARSEOps does not support arbitrary Ordinal types.  This is
    /// why this class' other_type typedef has a default definition
    /// which results in a compile error if used.  The graph<int,N>
    /// specialization (see below) is defined.
    template <class O, class N>
    struct graph {
      // This will raise a compile error if the given ordinal type O is not supported.
      typedef typename O::this_ordinal_not_supported_by_cusparse graph_type;
    };

    //! Partial specialization of graph<O,N> for O=int.
    template <class N>
    struct graph<int,N> {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    /// \brief Typedef for local matrix class.
    ///
    /// CUSPARSEOps does not support arbitrary Ordinal types.  This is
    /// why this class' other_type typedef has a default definition
    /// which results in a compile error if used.  The matrix<S,int,N>
    /// specialization (see below) is defined.
    template <class S, class O, class N>
    struct matrix {
      // This will raise a compile error if the given ordinal type O is not supported.
      typedef typename O::this_ordinal_not_supported_by_cusparse matrix_type;
    };

    //! Partial specialization of matrix<S,O,N> for O=int.
    template <class S, class N>
    struct matrix<S,int,N> {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar
    /// type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c
    /// CUSPARSEOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    CUSPARSEOps (const RCP<Node> &node) :
      Details::CUSPARSEOpsStub<Scalar, Node> (node) {}

    /// \brief "Sum constructor": compute *this = alpha*A + beta*B.
    ///
    /// The resulting matrix shares the Node instance and copies the
    /// parameters of the matrix A.
    CUSPARSEOps (const Scalar& alpha,
                 const CUSPARSEOps<Scalar, Node>& A,
                 const Scalar& beta,
                 const CUSPARSEOps<Scalar, Node>& B) :
      Details::CUSPARSEOpsStub<Scalar, Node> (alpha, A, beta, B) {}

    //! Destructor
    virtual ~CUSPARSEOps () {}

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosClassic::CUSPARSEOps<"
         << "Scalar=" << TypeNameTraits<Scalar>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //@}
  private:
    //! Copy constructor (protected and unimplemented)
    CUSPARSEOps (const CUSPARSEOps& source);
  };


  //! Partial specialization of CUSPARSEOps for Scalar=float.
  template <class Node>
  class CUSPARSEOps<float, Node> : public Details::CUSPARSEOpsGeneric<float, Node> {
  public:
    typedef float scalar_type;
    typedef int ordinal_type;
    typedef Node node_type;
    typedef CUSPARSEOps<float, Node> sparse_ops_type;

    template <class O, class N>
    struct graph {
      typedef typename O::this_ordinal_not_supported_by_cusparse graph_type;
    };

    template <class N>
    struct graph<int,N> {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    template <class S, class O, class N>
    struct matrix {
      typedef typename O::this_ordinal_not_supported_by_cusparse matrix_type;
    };

    template <class S, class N>
    struct matrix<S,int,N> {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    CUSPARSEOps (const RCP<Node> &node) :
      Details::CUSPARSEOpsGeneric<float, Node> (node) {}

    CUSPARSEOps (const float& alpha,
                 const CUSPARSEOps<float, Node>& A,
                 const float& beta,
                 const CUSPARSEOps<float, Node>& B) :
      Details::CUSPARSEOpsGeneric<float, Node> (alpha, A, beta, B) {}

    virtual ~CUSPARSEOps () {}

    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosClassic::CUSPARSEOps<"
         << "Scalar=float"
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

  private:
    CUSPARSEOps (const CUSPARSEOps& source);
  };


  //! Partial specialization of CUSPARSEOps for Scalar=double.
  template <class Node>
  class CUSPARSEOps<double, Node> : public Details::CUSPARSEOpsGeneric<double, Node> {
  public:
    typedef double scalar_type;
    typedef int ordinal_type;
    typedef Node node_type;
    typedef CUSPARSEOps<double, Node> sparse_ops_type;

    template <class O, class N>
    struct graph {
      typedef typename O::this_ordinal_not_supported_by_cusparse graph_type;
    };

    template <class N>
    struct graph<int,N> {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    template <class S, class O, class N>
    struct matrix {
      typedef typename O::this_ordinal_not_supported_by_cusparse matrix_type;
    };

    template <class S, class N>
    struct matrix<S,int,N> {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    CUSPARSEOps (const RCP<Node> &node) :
      Details::CUSPARSEOpsGeneric<double, Node> (node) {}

    CUSPARSEOps (const double& alpha,
                 const CUSPARSEOps<double, Node>& A,
                 const double& beta,
                 const CUSPARSEOps<double, Node>& B) :
      Details::CUSPARSEOpsGeneric<double, Node> (alpha, A, beta, B) {}

    virtual ~CUSPARSEOps () {}

    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosClassic::CUSPARSEOps<"
         << "Scalar=double"
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

  private:
    CUSPARSEOps (const CUSPARSEOps& source);
  };


  // Partial specialization for Scalar=void.  It omits instance
  // methods relating to Scalar, and methods relating to allocating a
  // matrix.
  template <class Node>
  class CUSPARSEOps<void, Node> : public Details::CUSPARSEOpsBase<Node> {
  public:
    //@{
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef void    scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef int     ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef CUSPARSEOps<void, Node> sparse_ops_type;

    /// \brief Typedef for local graph class.
    ///
    /// CUSPARSEOps does not support arbitrary Ordinal types.  This is
    /// why this class' other_type typedef has a default definition
    /// which results in a compile error if used.  The graph<int,N>
    /// specialization (see below) is defined.
    template <class O, class N>
    struct graph {
      // This will raise a compile error if the given ordinal type O is not supported.
      typedef typename O::this_ordinal_not_supported_by_cusparse graph_type;
    };

    //! Partial specialization of graph<O,N> for O=int.
    template <class N>
    struct graph<int,N> {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    /// \brief Typedef for local matrix class.
    ///
    /// CUSPARSEOps does not support arbitrary Ordinal types.  This is
    /// why this class' other_type typedef has a default definition
    /// which results in a compile error if used.  The matrix<S,int,N>
    /// specialization (see below) is defined.
    template <class S, class O, class N>
    struct matrix {
      // This will raise a compile error if the given ordinal type O is not supported.
      typedef typename O::this_ordinal_not_supported_by_cusparse matrix_type;
    };

    //! Partial specialization of matrix<S,O,N> for O=int.
    template <class S, class N>
    struct matrix<S,int,N> {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar
    /// type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c
    /// CUSPARSEOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    CUSPARSEOps (const RCP<Node> &node) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "CUSPARSEOps: Not allowed to instantiate this class for Scalar=void.");
    }

    //! Destructor
    virtual ~CUSPARSEOps () {}

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    virtual std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosClassic::CUSPARSEOps<"
         << "Scalar=void"
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    CUSPARSEOps (const CUSPARSEOps& source);
  };

} // namespace KokkosClassic

#endif /* KOKKOS_CUSPARSEOPS_HPP */

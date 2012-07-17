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

/*

notes:
* multiply requires encapsulation of multivector using make_array2d_view
* crs construction requires encapsulation of CSR matrix data via make_csr_matrix_view
  this object can be used for multiplication or conversion to another format.
* cusp::transpose() and cusp::convert() will come into play
*/

#ifndef KOKKOS_CUSPOPS_HPP
#define KOKKOS_CUSPOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_TypeTraits.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_CUDANodeUtils.hpp>

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"

#include <cusp/csr_matrix.h>
#include "Kokkos_CuspOps.cuh"

namespace Kokkos {


  //! \class CuspCrsGraph
  /** \brief CRS sparse graph class supporting the Cusp library.
  */
  template <class Ordinal, 
            class Node>
  class CuspCrsGraph : public CrsGraphBase<Ordinal,Node>
  {
    public:
      CuspCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params);
      bool isEmpty() const;
      void setStructure(const ArrayRCP<const Ordinal>  &ptrs,
                        const ArrayRCP<const Ordinal> &inds);
      void setDeviceData(const ArrayRCP<const Ordinal> &devptrs, const ArrayRCP<const Ordinal> &devinds);
      inline ArrayRCP<const Ordinal> getPointers() const;
      inline ArrayRCP<const Ordinal> getIndices() const;
      inline ArrayRCP<const Ordinal> getDevPointers() const;
      inline ArrayRCP<const Ordinal> getDevIndices() const;
      inline bool isInitialized() const;
      void setMatDesc(Teuchos::EUplo  uplo, Teuchos::EDiag  diag);
      void getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;
    private:
      bool isInitialized_;
      bool isEmpty_;
      Teuchos::EUplo uplo_;
      Teuchos::EDiag diag_;
      // graph data
      ArrayRCP<const Ordinal> host_rowptrs_, dev_rowptrs_, host_colinds_, dev_colinds_;
  };

  //! \class CuspCrsMatrix
  /** \brief CRS sparse matrix class supporting the Cusp library.
  */
  template <class Scalar,
            class Ordinal,
            class Node>
  class CuspCrsMatrix : public CrsMatrixBase<Scalar,Ordinal,Node>
  {
    public:
      CuspCrsMatrix(const RCP<const CuspCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params);
      void setValues(const ArrayRCP<const Scalar> &vals);
      void setDeviceData(const ArrayRCP<const Scalar> &devvals);
      inline ArrayRCP<const Scalar> getValues() const;
      inline ArrayRCP<const Scalar> getDevValues() const;
      inline bool isInitialized() const;
    private:
      bool isInitialized_;
      // matrix data
      ArrayRCP<const Scalar> vals_, dev_vals_;
      // TODO: add CSC data, for efficient transpose multiply
  };

  template <class Ordinal, class Node>
  CuspCrsGraph<Ordinal,Node>::CuspCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params)
  : CrsGraphBase<Ordinal,Node>(numRows,numCols,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Ordinal, class Node>
  bool CuspCrsGraph<Ordinal,Node>::isEmpty() const
  { return isEmpty_; }

  template <class Ordinal, class Node>
  void CuspCrsGraph<Ordinal,Node>::setStructure(const ArrayRCP<const Ordinal> &ptrs,
                                        const ArrayRCP<const Ordinal> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const Ordinal numrows = this->getNumRows();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)ptrs.size() != (size_t)numrows+1
        || ptrs[0] != 0
        || (size_t)inds.size() != (size_t)ptrs[numrows],
        std::runtime_error, " graph data not coherent."
    )
    const Ordinal numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix has already been initialized"
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    host_rowptrs_ = ptrs;
    host_colinds_ = inds;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> CuspCrsGraph<Ordinal,Node>::getPointers() const
  { return host_rowptrs_; }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> CuspCrsGraph<Ordinal,Node>::getIndices() const
  { return host_colinds_; }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> CuspCrsGraph<Ordinal,Node>::getDevPointers() const
  { return dev_rowptrs_; }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> CuspCrsGraph<Ordinal,Node>::getDevIndices() const
  { return dev_colinds_; }

  template <class Ordinal, class Node>
  void CuspCrsGraph<Ordinal,Node>::setDeviceData(const ArrayRCP<const Ordinal> &devptrs,
                                         const ArrayRCP<const Ordinal> &devinds)
  { dev_rowptrs_ = devptrs; dev_colinds_ = devinds; }

  template <class Ordinal, class Node>
  bool CuspCrsGraph<Ordinal,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Ordinal, class Node>
  void CuspCrsGraph<Ordinal,Node>::setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    std::string tfecfFuncName("setMatDesc()");
    uplo_ = uplo;
    diag_ = diag;
  }

  template <class Ordinal, class Node>
  void CuspCrsGraph<Ordinal,Node>::getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = uplo_;
    diag = diag_;
  }

  template <class Scalar, class Ordinal, class Node>
  CuspCrsMatrix<Scalar,Ordinal,Node>::CuspCrsMatrix(const RCP<const CuspCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,Ordinal,Node>(graph,params)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template <class Scalar, class Ordinal, class Node>
  void CuspCrsMatrix<Scalar,Ordinal,Node>::setValues(const ArrayRCP<const Scalar> &vals)
  {
    std::string tfecfFuncName("setValues(vals)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix is already initialized."
    )
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<const Scalar> CuspCrsMatrix<Scalar,Ordinal,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Ordinal, class Node>
  bool CuspCrsMatrix<Scalar,Ordinal,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<const Scalar> CuspCrsMatrix<Scalar,Ordinal,Node>::getDevValues() const
  { return dev_vals_; }

  template <class Scalar, class Ordinal, class Node>
  void CuspCrsMatrix<Scalar,Ordinal,Node>::setDeviceData(const ArrayRCP<const Scalar> &devvals)
  { dev_vals_ = devvals; }

  /// \class CuspOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Ordinal, class Node>
  class CuspOps {
  public:
    //@{
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef CuspOps<Scalar,Ordinal,Node> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef CuspCrsGraph<O,N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef CuspCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c
    /// CuspOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CuspOps<S2,Ordinal,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    CuspOps(const RCP<Node> &node);

    //! Destructor
    ~CuspOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &rowPtrs);

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const Ordinal> &ptrs);

    //! Finalize a graph is null for Cusp.
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, CuspCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const CuspCrsGraph<Ordinal,Node> &graph, CuspCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, CuspCrsGraph<Ordinal,Node> &graph, CuspCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    void setGraphAndMatrix(const RCP<const CuspCrsGraph<Ordinal,Node> > &graph,
                           const RCP<const CuspCrsMatrix<Scalar,Ordinal,Node> > &matrix);

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
              MultiVector<RangeScalar,Node> &Y) const;

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
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result: \f$Y := Y + \alpha A X\f$.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

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
           MultiVector<RangeScalar,Node> &X) const {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Kokkos::CuspOps does not support solve() because Cusp library does not provide sparse triangular solve.");
    }

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    CuspOps(const CuspOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    Ordinal numRows_, numCols_, numNZ_;
    bool isInitialized_;

    typedef typename Cuspdetails::types<Ordinal,Scalar>::csr_view           csr_matrix_view;

    ArrayRCP<const Ordinal> rowPtrs_, colInds_;
    ArrayRCP<const Scalar>  rowVals_;
    RCP<csr_matrix_view>    A_mat_, transA_mat;
  };


  // ======= matrix finalization ===========
  template <class Scalar, class Ordinal, class Node>
  void CuspOps<Scalar,Ordinal,Node>::finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag,
                                                   CuspCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params)
  {
    const std::string prefix("finalizeGraph()");
    RCP<Node> node = graph.getNode();
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, prefix << ": graph has not yet been initialized."
    )
    // diag: have to allocate and indicate
    ArrayRCP<Ordinal>       devinds, devptrs;
    ArrayRCP<const Ordinal> hostinds = graph.getIndices();
    ArrayRCP<const Ordinal> hostptrs = graph.getPointers();
    const Ordinal numRows = graph.getNumRows();
    // set description
    graph.setMatDesc(uplo,diag);
    // allocate and initialize data
    devptrs = node->template allocBuffer<Ordinal>( numRows+1 );
    node->copyToBuffer( numRows+1,  hostptrs(),  devptrs );
    const Ordinal numnz = hostinds.size();
    if (numnz) {
      devinds = node->template allocBuffer<Ordinal>( numnz );
      node->copyToBuffer( numnz,     hostinds(), devinds );
    }
    // set the data
    graph.setDeviceData(devptrs,devinds);
  }

  // ======= matrix finalization ===========
  template <class Scalar, class Ordinal, class Node>
  void CuspOps<Scalar,Ordinal,Node>::
  finalizeMatrix (const CuspCrsGraph<Ordinal,Node> &graph,
                        CuspCrsMatrix<Scalar,Ordinal,Node> &matrix,
                  const RCP<ParameterList> &params)
  {
    std::string FuncName("Kokkos::CuspOps::finalizeMatrix()");
    RCP<Node> node = graph.getNode();
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
    // diag: have to allocate and indicate
    Teuchos::EUplo uplo;
    Teuchos::EDiag diag;
    graph.getMatDesc(uplo,diag);
    ArrayRCP<Scalar>           devvals;
    ArrayRCP<const Ordinal>    hostptrs = graph.getPointers();
    ArrayRCP<const Scalar> hostvals = matrix.getValues();
    const Ordinal numRows = graph.getNumRows();
    const Ordinal numnz = hostptrs[numRows];
    if (numnz) {
      devvals = node->template allocBuffer<Scalar>( numnz );
      node->copyToBuffer( numnz,     hostvals(), devvals );
    }
    matrix.setDeviceData(devvals);
    if (params != null) {
      // warn? if (params->get("Prepare Solve",false))  {}
      // warn? if (params->get("Prepare Tranpose Solve",false))  {}
      // warn? if (params->get("Prepare Conjugate Tranpose Solve",false))  {}
      if (params->get("Prepare Transpose Multiply",false))  {} // explicitly create tranpose if true
      if (params->get("Prepare Multiply",true))  {}            // delete non-transpose if false
    }
  }

  // ======= graph and matrix finalization ===========
  template <class Scalar, class Ordinal, class Node>
  void CuspOps<Scalar,Ordinal,Node>::finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag,
                                                            CuspCrsGraph<Ordinal,Node> &graph, 
                                                            CuspCrsMatrix<Scalar,Ordinal,Node> &matrix,
                                                            const RCP<ParameterList> &params)
  {
    std::string FuncName("Kokkos::CuspOps::finalizeGraphAndMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
    // no benefit to doing them together; do them separately
    finalizeGraph(uplo,diag,graph,params);
    finalizeMatrix(graph,matrix,params);
  }


  // ======= pointer allocation ===========
  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<Ordinal>
  CuspOps<Scalar,Ordinal,Node>::allocRowPtrs(const RCP<Node> &/*node*/, const ArrayView<const Ordinal> &numEntriesPerRow)
  {
    // alloc page-locked ("pinned") memory on the host, specially allocated and specially deallocated
    CUDANodeHostPinnedDeallocator<Ordinal> dealloc;
    ArrayRCP<Ordinal> ptrs = dealloc.alloc(numEntriesPerRow.size() + 1);
    ptrs[0] = 0;
    std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
    return ptrs;
  }

  // ======= other allocation ===========
  template <class Scalar, class Ordinal, class Node>
  template <class T>
  ArrayRCP<T>
  CuspOps<Scalar,Ordinal,Node>::allocStorage(const RCP<Node> &/*node*/, const ArrayView<const Ordinal> &rowPtrs)
  {
    // alloc page-locked ("pinned") memory on the host, specially allocated and specially deallocated
    const Ordinal totalNumEntries = *(rowPtrs.end()-1);
    CUDANodeHostPinnedDeallocator<T> dealloc;
    ArrayRCP<T> buf = dealloc.alloc(totalNumEntries);
    std::fill(buf.begin(), buf.end(), Teuchos::ScalarTraits<T>::zero() );
    return buf;
  }

  template<class Scalar, class Ordinal, class Node>
  CuspOps<Scalar,Ordinal,Node>::CuspOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , numCols_(0)
  , numNZ_(0)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are CUDA Nodes
    Teuchos::CompileTimeAssert<Node::isCUDANode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node>
  CuspOps<Scalar,Ordinal,Node>::~CuspOps()
  { }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> CuspOps<Scalar,Ordinal,Node>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node>
  void CuspOps<Scalar,Ordinal,Node>::setGraphAndMatrix(const RCP<const CuspCrsGraph<Ordinal,Node> > &graph, 
                                                       const RCP<const CuspCrsMatrix<Scalar,Ordinal,Node> > &matrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " operators already initialized.");
    // get data from the matrix
    numRows_ = graph->getNumRows();
    numRows_ = graph->getNumCols();
    rowPtrs_ = graph->getDevPointers();
    colInds_ = graph->getDevIndices();
    rowVals_ = matrix->getDevValues();
    numNZ_ = colInds_.size();
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspOps<Scalar,Ordinal,Node>::solve(Teuchos::ETransp trans,
                                           const MultiVector<DomainScalar,Node> &Y,
                                                 MultiVector< RangeScalar,Node> &X) const
  {
    //
    std::string tfecfFuncName("solve()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    Ordinal stride_x = (Ordinal)X.getStride(),
            stride_y = (Ordinal)Y.getStride();
    const Scalar * data_y = Y.getValues().getRawPtr();
    Scalar * data_x = X.getValuesNonConst().getRawPtr();
    const Ordinal numMatRows = X.getNumRows();
    const Ordinal numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumRows() != Y.getNumRows(),
        std::runtime_error, "X and Y do not have the same number of row vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        Y.getNumRows() != numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspOps<Scalar,Ordinal,Node>::multiply(Teuchos::ETransp trans,
                                          RangeScalar alpha,
                                          const MultiVector<DomainScalar,Node> &X,
                                                MultiVector< RangeScalar,Node> &Y) const
  {
    //
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    Ordinal stride_x = (Ordinal)X.getStride(),
            stride_y = (Ordinal)Y.getStride();
    const Scalar * data_x = X.getValues().getRawPtr();
    Scalar * data_y = Y.getValuesNonConst().getRawPtr();
    const Ordinal numMatRows = Y.getNumRows();
    const Ordinal numMatCols = X.getNumRows();
    const Ordinal numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)Y.getNumRows() != (size_t)numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspOps<Scalar,Ordinal,Node>::multiply(Teuchos::ETransp trans,
                                          RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                          RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    //
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, "sparse operators have not been initialized with graph and matrix data; call setGraphAndMatrix() first.")
    // get pointers,stride from X and Y
    Ordinal stride_x = (Ordinal)X.getStride(),
            stride_y = (Ordinal)Y.getStride();
    const Scalar * data_x = X.getValues().getRawPtr();
    Scalar * data_y = Y.getValuesNonConst().getRawPtr();
    const Ordinal numMatRows = Y.getNumRows();
    const Ordinal numMatCols = X.getNumRows();
    const Ordinal numRHS     = X.getNumCols();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, "X and Y do not have the same number of column vectors.")
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)Y.getNumRows() != (size_t)numRows_,
        std::runtime_error, "Y does not have the same number of rows as does the matrix.")
  }

} // namespace Kokkos

#endif /* KOKKOS_CUSPOPS_HPP */

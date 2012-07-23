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

#ifndef KOKKOS_DEFAULTSPARSEOPS_HPP
#define KOKKOS_DEFAULTSPARSEOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseSolveKernelOps.hpp"
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

namespace Kokkos {

  namespace details {

    template <class O>
    struct rowPtrsInitKernel {
      O *rowPtrs;
      inline KERNEL_PREFIX void execute(int i) const {
        rowPtrs[i] = 0;
      }
    };

    template <class O, class T>
    struct valsInitKernel {
      const O *rowPtrs;
      T * entries;
      inline KERNEL_PREFIX void execute(int i) const {
        T *beg = entries+rowPtrs[i],
          *end = entries+rowPtrs[i+1];
        while (beg != end) {
          *beg++ = Teuchos::ScalarTraits<T>::zero();
        }
      }
    };

    template <class Ordinal, class Node>
    class FirstTouchCRSAllocator {
      public:
      //! \brief Allocate and initialize the storage for the matrix values.
      static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &numEntriesPerRow)
      {
        const Ordinal numrows = numEntriesPerRow.size();
        details::rowPtrsInitKernel<Ordinal> kern;
        // allocate
        kern.rowPtrs = new Ordinal[numrows+1];
        // parallel first touch
        node->parallel_for(0,numrows+1,kern);
        // encapsulate
        ArrayRCP<Ordinal> ptrs = arcp<Ordinal>(kern.rowPtrs,0,numrows+1,true);
        // compute in serial. parallelize later, perhaps; it's only O(N)
        ptrs[0] = 0;
        std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
        return ptrs;
      }

      //! \brief Allocate and initialize the storage for a sparse graph.
      template <class T>
      static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const Ordinal> &rowPtrs)
      {
        const Ordinal totalNumEntries = *(rowPtrs.end()-1);
        const Ordinal numRows = rowPtrs.size() - 1;
        details::valsInitKernel<Ordinal,T> kern;
        ArrayRCP<T> vals;
        if (totalNumEntries > 0) {
          // allocate
          kern.entries = new T[totalNumEntries];
          kern.rowPtrs = rowPtrs.getRawPtr();
          // first touch
          node->parallel_for(0,numRows,kern);
          // encapsulate
          vals = arcp<T>(kern.entries,0,totalNumEntries,true);
        }
        return vals;
      }
    };

    template <class Ordinal, class Node>
    class DefaultCRSAllocator {
      public:
      //! \brief Allocate and initialize the storage for the matrix values.
      static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &numEntriesPerRow)
      {
        ArrayRCP<Ordinal> ptrs = arcp<Ordinal>( numEntriesPerRow.size() + 1 );
        ptrs[0] = 0;
        std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
        return ptrs;
      }

      //! \brief Allocate and initialize the storage for a sparse graph.
      template <class T>
      static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const Ordinal> &rowPtrs)
      {
        const Ordinal totalNumEntries = *(rowPtrs.end()-1);
        // alloc data
        ArrayRCP<T> vals;
        if (totalNumEntries > 0) vals = arcp<T>(totalNumEntries);
        std::fill( vals.begin(), vals.end(), Teuchos::ScalarTraits<T>::zero() );
        return vals;
      }
    };

  }

  //! \class DefaultCrsGraph
  /** \brief Default implementation of CRS sparse graph, using generic kernels and suitable for host-based nodes.
  */
  template <class Ordinal,
            class Node>
  class DefaultCrsGraph : public CrsGraphBase<Ordinal,Node>
  {
    public:
      DefaultCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params);
      bool isEmpty() const;
      void setStructure(const ArrayRCP<const Ordinal> &ptrs,
                        const ArrayRCP<const Ordinal> &inds);
      inline ArrayRCP<const Ordinal> getPointers() const;
      inline ArrayRCP<const Ordinal> getIndices() const;
      inline bool isInitialized() const;
      inline void setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag);
      inline void getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;
    private:
      ArrayRCP<const Ordinal> ptrs_;
      ArrayRCP<const Ordinal> inds_;
      bool isInitialized_;
      bool isEmpty_;
      Teuchos::EUplo  tri_uplo_;
      Teuchos::EDiag unit_diag_;
  };

  //! \class DefaultCrsMatrix
  /** \brief Default implementation of CRS sparse matrix, using generic kernels and suitable for host-based nodes.
  */
  template <class Scalar,
            class Ordinal,
            class Node>
  class DefaultCrsMatrix : public CrsMatrixBase<Scalar,Ordinal,Node>
  {
    public:
      DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params);
      void setValues(const ArrayRCP<const Scalar> &vals);
      inline ArrayRCP<const Scalar> getValues() const;
      inline bool isInitialized() const;
    private:
      ArrayRCP<const Scalar> vals_;
      bool isInitialized_;
  };

  template <class Ordinal, class Node>
  DefaultCrsGraph<Ordinal,Node>::DefaultCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params)
  : CrsGraphBase<Ordinal,Node>(numRows,numCols,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isEmpty() const
  { return isEmpty_; }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setStructure(
                      const ArrayRCP<const Ordinal> &ptrs,
                      const ArrayRCP<const Ordinal> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const Ordinal numrows = this->getNumRows();

    // mfh 19 June 2012: The tests expect std::runtime_error rather
    // than the arguably more appropriate std::invalid_argument, so
    // I'll throw std::runtime_error here.  Ditto for the other checks
    // below.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptrs.is_null (),
      std::runtime_error,
      ": The input array 'ptrs' must be nonnull, even for a matrix with zero "
      "rows.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptrs.size() != (size_t) numrows+1,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptrs.size() = " << ptrs.size() << " != numrows+1 = "
      << (numrows+1) << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptrs[0] != 0,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptrs[0] = " << ptrs[0] << " != 0.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) inds.size() != (size_t) ptrs[numrows],
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- inds.size() = " << inds.size() << " != ptrs[numrows="
      << numrows << "] = " << ptrs[numrows] << ".");

    const Ordinal numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_,
      std::runtime_error,
      " matrix has already been initialized."
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    ptrs_ = ptrs;
    inds_ = inds;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsGraph<Ordinal,Node>::getPointers() const
  { return ptrs_; }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsGraph<Ordinal,Node>::getIndices() const
  { return inds_; }


  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    tri_uplo_ = uplo;
    unit_diag_ = diag;
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = tri_uplo_;
    diag = unit_diag_;
  }

  template <class Scalar, class Ordinal, class Node>
  DefaultCrsMatrix<Scalar,Ordinal,Node>::DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,Ordinal,Node>(graph,params)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultCrsMatrix<Scalar,Ordinal,Node>::setValues(const ArrayRCP<const Scalar> &vals)
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
  ArrayRCP<const Scalar> DefaultCrsMatrix<Scalar,Ordinal,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Ordinal, class Node>
  bool DefaultCrsMatrix<Scalar,Ordinal,Node>::isInitialized() const
  {
    return isInitialized_;
  }

  /// \class DefaultHostSparseOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Ordinal, class Node, class Allocator = details::DefaultCRSAllocator<Ordinal,Node> >
  class DefaultHostSparseOps {
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
    typedef DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
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
    /// DefaultHostSparseOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef DefaultHostSparseOps<S2,Ordinal,Node,Allocator> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    DefaultHostSparseOps(const RCP<Node> &node);

    //! Constructor accepting and retaining a node object, and taking parameters.
    DefaultHostSparseOps(const RCP<Node> &node, Teuchos::ParameterList& params);

    //! Destructor
    ~DefaultHostSparseOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &numEntriesPerRow)
    {
      return Allocator::allocRowPtrs(node,numEntriesPerRow);
    }

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const Ordinal> &rowPtrs)
    {
      return Allocator::template allocStorage<T>(node,rowPtrs);
    }

    //! Finalize a graph
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    ///
    /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
    ///   else LOWER_TRI if the matrix is lower triangular.
    ///
    /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
    ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
    ///
    void setGraphAndMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph,
                           const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &matrix);

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
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector. Contents will be overwritten.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := beta * Y + alpha * Op(A) * X.
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
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param beta [in] Scalar constant \f$\beta\f$ by which to
    ///   multiply Y when summing with the result of the sparse
    ///   matrix-(multi)vector multiply.
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
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultHostSparseOps(const DefaultHostSparseOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    // we do this one of two ways:
    // packed CRS: array of row pointers, array of indices, array of values.
    ArrayRCP<const Ordinal> inds_;
    ArrayRCP<const Ordinal> ptrs_;
    ArrayRCP<const Scalar>  vals_;

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;

    Ordinal numRows_;
    bool isInitialized_;
    bool isEmpty_;
  };


  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params)
  {
    // nothing much to do here
    graph.setMatDesc(uplo,diag);
    std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeGraph(graph,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // nothing much to do here
    std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // nothing much to do here
    graph.setMatDesc(uplo,diag);
    std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeGraphAndMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }


  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::DefaultHostSparseOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::
  DefaultHostSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params)
  : node_(node)
  , numRows_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    (void) params; // Not using this yet.

    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::~DefaultHostSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  RCP<Node> DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::setGraphAndMatrix(
              const RCP<const DefaultCrsGraph<Ordinal,Node> > &opgraph,
              const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " operators already initialized.");
    numRows_ = opgraph->getNumRows();
    if (opgraph->isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      ptrs_ = opgraph->getPointers();
      inds_ = opgraph->getIndices();
      vals_ = opmatrix->getValues();
      // these checks just about the most that we can perform
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)ptrs_.size() != (size_t)numRows_+1
          || inds_.size() != vals_.size(),
          std::runtime_error, " matrix and graph seem incongruent.");
    }
    opgraph->getMatDesc( tri_uplo_, unit_diag_ );
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::solve(Teuchos::ETransp trans,
                                                        const MultiVector<DomainScalar,Node> &Y,
                                                              MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    typedef DefaultSparseSolveOp<Scalar,Ordinal,DomainScalar,RangeScalar>            Op;
    typedef DefaultSparseTransposeSolveOp<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " this solve was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumCols() != (size_t)Y.getNumCols(),
        std::runtime_error, " Left hand side and right hand side multivectors have differing numbers of vectors."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumRows() < (size_t)numRows_,
        std::runtime_error, " Left-hand-side multivector does not have enough rows. "
                            "Likely cause is that the column map was not provided to "
                            "the Tpetra::CrsMatrix in the case of an implicit unit diagonal."
    );

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(unit_diag_ != Teuchos::UNIT_DIAG, std::runtime_error,
          " solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<    Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<    Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X,
                                      MultiVector<RangeScalar ,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    // the 1 template parameter below means that beta is not used in computations
    // and the output multivector enjoys overwrite semantics (i.e., will overwrite data/NaNs in Y)
    typedef DefaultSparseMultiplyOp<         Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op;
    typedef DefaultSparseTransposeMultiplyOp<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<     Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<     Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    // the 0 template parameter below means that beta is used in computations
    // and the output multivector enjoys accumulation semantics (i.e., will not overwrite data/NaNs in Y)
    typedef DefaultSparseMultiplyOp<         Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op;
    typedef DefaultSparseTransposeMultiplyOp<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<     Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.ptrs    = rbh.template addConstBuffer<     Ordinal>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }


  /** \example CrsMatrix_DefaultMultiplyTests.hpp
    * This is an example that unit tests and demonstrates the implementation requirements for the DefaultSparseOps class.
    */

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */


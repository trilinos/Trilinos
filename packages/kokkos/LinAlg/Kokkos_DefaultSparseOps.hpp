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

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseScaleKernelOps.hpp"
#include "Kokkos_DefaultSparseSolveKernelOps.hpp"
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

namespace Kokkos {

  // forward declaration of DefaultCrsMatrix
  template <class Scalar, class Ordinal, class Node> class DefaultCrsMatrix;

  //! \class DefaultCrsGraph 
  /** \brief Default implementation of CRS sparse graph, using generic kernels and suitable for host-based nodes.
  */
  template <class Ordinal, 
            class Node>
  class DefaultCrsGraph : public CrsGraphBase<Ordinal,Node> 
  {
    public:
      DefaultCrsGraphBase(size_t numRows, const RCP<Node> &node);
      virtual size_t getNumEntries() const;
      virtual bool isEmpty() const;
      virtual bool isFinalized() const;
      virtual void setStructure(const ArrayRCP<const size_t>  &ptrs,
                                const ArrayRCP<const Ordinal> &inds);
      void finalize(Teuchos::ParamterList &params);
      void finalizeMatrix(DefaultCrsMatrix<Scalar,Ordinal,Node> &mat, Teuchos::ParamterList &params) const;
      void finalizeGraphAndMatrix(DefaultCrsMatrix<Scalar,Ordinal,Node> &mat, Teuchos::ParamterList &params);
      inline ArrayRCP<const size_t> getPointers() const;
      inline ArrayRCP<const Ordinal> getIndices() const;
    private:
      ArrayRCP<const size_t>  ptrs_;
      ArrayRCP<const Ordinal> inds_;
      size_t numEntries_;
      bool isInitialized_, isEmpty_, isFinalized_;
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
      DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph);
      virtual void setValues(const ArrayRCP<const Scalar> &vals);
      inline ArrayRCP<const Scalar> getValues() const;
      bool isInitialized() const;
    private:
      ArrayRPC<const Scalar> values_;
      bool isInitialized_;
  };

  template <class Ordinal, class Node>
  DefaultCrsGraph<Ordinal,Node>::DefaultCrsGraphBase(size_t numRows, const RCP<Node> &node)
  : CrsGraphBase<Ordinal,Node>(numRows,node)
  , numEntries_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  , isFinalized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Ordinal, class Node>
  size_t DefaultCrsGraph<Ordinal,Node>::getNumEntries() const
  {
    return numEntries_;
  }

  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isEmpty() const
  {
    return isEmpty_;
  }

  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isFinalized() const
  {
    return isFinalized_;
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setStructure(
                      const ArrayRCP<const size_t>  &ptrs,
                      const ArrayRCP<const Ordinal> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ptrs.size() != getNumRows()+1 
        || ptrs[0] != 0
        || inds.size() != ptrs[getNumRows()],
        std::runtime_error, " graph data not coherent."
    )
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix has already been initialized"
    )
    if (numRows_ == 0 || numEntries_ == 0) isEmpty_ = true;
    numEntries_ = ptrs[getNumRows()];
    ptrs_ = ptrs;
    inds_ = inds;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::finalize(Teuchos::ParamterList &params)
  {
    finalized_ = true;
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::finalizeMatrix(
                      DefaultCrsMatrix<Scalar,Ordinal,Node> &mat, 
                      Teuchos::ParamterList &params) const
  {
    // not much to do here
    std::string tfecfFuncName("finalizeMatrix(matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        mat.isInitialized() == false,
        std::runtime_error, " matrix not initialized yet."
    )
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::finalizeGraphAndMatrix(
                      DefaultCrsMatrix<Scalar,Ordinal,Node> &mat, 
                      Teuchos::ParamterList &params) 
  { 
    // not much to do here
    std::string tfecfFuncName("finalizeGraphAndMatrix(matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " graph not initialized yet."
    )
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        mat.isInitialized() == false,
        std::runtime_error, " matrix not initialized yet."
    )
    isFinalized_ = true;
  }

  template <class Ordinal, class Node>
  ArrayRCP<const size_t> DefaultCrsGraph<Ordinal,Node>::getPointers() const
  {
    return ptrs_:
  }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsGraph<Ordinal,Node>::getIndices() const
  {
    return inds_;
  }

  template <class Scalar, class Ordinal, class Node>
  DefaultCrsMatrix<Scalar,Ordinal,Node>::DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph)
  : CrsMatrixBase<Scalar,Ordinal,Node>(graph) 
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
        vals.size() != graph_->getNumEntries(),      
        std::runtime_error, " provided values are not congruent with graph structure."
    )
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsMatrix<Scalar,Ordinal,Node>::getIndices() const
  { return vals_; }

  /// \class DefaultHostSparseOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Ordinal, class Node>
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
    typedef DefaultHostSparseOps<Scalar,Ordinal,Node> sparse_ops_type;

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
    /// The "rebind" struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// This always specifies a specialization of \c
    /// DefaultHostSparseOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct rebind {
      typedef DefaultHostSparseOps<S2,Ordinal,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    DefaultHostSparseOps(const RCP<Node> &node);

    //! Destructor
    ~DefaultHostSparseOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Initialize sparse operations with a graph and matrix
    void setGraphAndMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, const DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix);

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
    /// \param diag [in] UNIT_DIAG if the matrix has unit diagonal,
    ///   else NON_UNIT_DIAG.
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           Teuchos::EUplo uplo,
           Teuchos::EDiag diag,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    /* Commented these out for now; the ability to change the data 
       means that the values can't be stored const; need to wait and see
       what the impact of that is */
    // //! Left-scale the matrix by the given vector.
    // template <class VectorScalar>
    // void leftScale(const MultiVector<VectorScalar,Node> &X);

    // //! Right-scale the matrix by the given vector.
    // template <class VectorScalar>
    // void rightScale(const MultiVector<VectorScalar,Node> &X);

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultHostSparseOps(const DefaultHostSparseOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    // we do this one of two ways:
    // packed CRS: array of row pointers, array of indices, array of values.
    ArrayRCP<const Ordinal> inds_;
    ArrayRCP<const size_t>  ptrs_;
    ArrayRCP<Scalar>  vals_;

    size_t numRows_;
    bool isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultHostSparseOps<Scalar,Ordinal,Node>::DefaultHostSparseOps(const RCP<Node> &node)
  : node_(node)
  {
    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultHostSparseOps<Scalar,Ordinal,Node>::~DefaultHostSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> DefaultHostSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::clear() {
    ptrs_          = null;
    inds_          = null;
    vals_          = null;
    numRows_       = 0;
    isInitialized_ = false;
    isEmpty_       = false;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::setGraphAndMatrix(const CrsGraph<Ordinal,Node> &graph, const CrsMatrix<Scalar,Ordinal,Node> &matrix);
  {
    std::string tfecfFuncName("setGraphAndMatrix(graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true, 
        std::runtime_error, " operators already initialized.");
    // can't do much more than this
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
        graph.isEmpty() == matrix.isEmpty() && graph.getNumRow() == matrix.getNumRows(), 
        std::runtime_error, " matrix and graph are not congruent.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      ptrs_ = graph.getPointers();
      inds_ = graph.getIndices();
      vals_ = matrix.getValues();
      // these checks just about the most that we can perform
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( 
          ptrs_.size() != numRows_+1 
          || inds_.size() != vals_.size(),
          std::runtime_error, " matrix and graph seem incongruent.");
    }
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag,
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,uplo,diag,Y,X)");
    typedef DefaultSparsePackedSolve<Scalar,Ordinal,DomainScalar,RangeScalar>            Op;
    typedef DefaultSparsePackedTransposeSolve<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(indsInit_ == false || valsInit_ == false, 
        std::runtime_error, " this solve was not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumCols() != Y.getNumCols(), 
        std::runtime_error, " Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumRows() < numRows_, 
        std::runtime_error, " Left-hand-side multivector does not have enough rows. "
                            "Likely cause is that the column map was not provided to "
                            "the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(diag != Teuchos::UNIT_DIAG, std::runtime_error,
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
        wdp.begs    = rbh.template addConstBuffer<     size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op>(0,numRHS,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.begs    = rbh.template addConstBuffer<     size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                      MultiVector<RangeScalar ,Node> &Y) const 
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    // the 1 template parameter below means that beta is not used in computations
    // and the output multivector enjoys overwrite semantics (i.e., will overwrite data/NaNs in Y)
    typedef DefaultSparsePackedMult<         Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op;
    typedef DefaultSparsePackedTransposeMult<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        indsInit_ == false || valsInit_ == false, 
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
        wdp.begs    = rbh.template addConstBuffer<      size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op>(0,numRows_*numRHS,wdp);
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
        wdp.begs    = rbh.template addConstBuffer<      size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    // the 0 template parameter below means that beta is used in computations
    // and the output multivector enjoys accumulation semantics (i.e., will not overwrite data/NaNs in Y)
    typedef DefaultSparsePackMult<     Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op;
    typedef DefaultSparsePackTransMult<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        indsInit_ == false || valsInit_ == false, 
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
        wdp.begs    = rbh.template addConstBuffer<      size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op>(0,numRows_*numRHS,wdp);
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
        wdp.begs    = rbh.template addConstBuffer<      size_t>(ptrs_);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp>(0,numRHS,wdp);
      }
    }
    return;
  }

//   template <class Scalar, class Ordinal, class Node>
//   template <class VectorScalar>
//   void DefaultHostSparseOps<Scalar,Ordinal,Node>::leftScale(const MultiVector<VectorScalar,Node> &X)
//   {
//     typedef DefaultSparseScaleOp1<Scalar,Ordinal,VectorScalar>  Op1D;
//     typedef DefaultSparseScaleOp2<Scalar,Ordinal,VectorScalar>  Op2D;
//     TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
//         Teuchos::typeName(*this) << "::leftScale(): operation not fully initialized.");
//     TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != 1);
//     ReadyBufferHelper<Node> rbh(node_);
//     if (begs1D_ != null) {
//       Op1D wdp;
//       rbh.begin();
//       wdp.numRows = numRows_;
//       wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
//       wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
//       wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
//       wdp.vals    = rbh.template addNonConstBuffer<Scalar>(vals1D_);
//       wdp.x       = rbh.template addConstBuffer<VectorScalar>(X.getValues());
//       rbh.end();
//       node_->template parallel_for<Op1D>(0,numRows_,wdp);
//     }
//     else {
//       Op2D wdp;
//       rbh.begin();
//       wdp.numRows = numRows_;
//       wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
//       wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
//       wdp.vals_beg   = rbh.template addNonConstBuffer<Scalar *>(valPtrs_);
//       wdp.x          = rbh.template addConstBuffer<VectorScalar>(X.getValues());
//       rbh.end();
//       node_->template parallel_for<Op2D>(0,numRows_,wdp);
//     }
//     return;
//   }

  /** \example CrsMatrix_DefaultMultiplyTests.hpp
    * This is an example that unit tests and demonstrates the implementation requirements for the DefaultSparseOps class.
    */

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */


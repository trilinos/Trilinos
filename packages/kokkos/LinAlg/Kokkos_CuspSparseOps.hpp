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

#ifndef KOKKOS_CUSPSPARSEOPS_HPP
#define KOKKOS_CUSPSPARSEOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp" 
#include "Kokkos_CrsGraphBase.hpp" 

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <CUSP/hyb_matrix.h>
#include <CUSP/csr_matrix.h>

namespace Kokkos {

  namespace CUSP {

    //! \class CUSPSparseGraph 
    /** \brief Container for CRS sparse graph for use with CUSPSparseOps.
    */
    template <class Ordinal, 
              class Node>
    class CUSPSparseGraph : public CrsGraphBase<Ordinal,Node> 
    {
      public:
        CUSPSparseGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params);
        bool isEmpty() const;
        void setStructure(const ArrayRCP<const size_t>  &ptrs,
                          const ArrayRCP<const Ordinal> &inds);
        inline ArrayRCP<const size_t> getPointers() const;
        inline ArrayRCP<const Ordinal> getIndices() const;
        inline bool isInitialized() const;
      private:
        ArrayRCP<const size_t>  ptrs_;
        ArrayRCP<const Ordinal> inds_;
        bool isInitialized_, isEmpty_;
    };
    
    //! \class CUSPSparseMatrix 
    /** \brief Container for sparse matrix, for use with CUSPSparseOps.
    */
    template <class Scalar, 
              class Ordinal, 
              class Node> 
    class CUSPSparseMatrix : public CrsMatrixBase<Scalar,Ordinal,Node> 
    {
      public:
        CUSPSparseMatrix(const RCP<const CUSPSparseGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params);
        void setValues(const ArrayRCP<const Scalar> &vals);
        inline ArrayRCP<const Scalar> getValues() const;
        inline bool isInitialized() const;
      private:
        ArrayRCP<const Scalar> vals_;
        bool isInitialized_;
    };
    
    template <class Ordinal, class Node>
    CUSPSparseGraph<Ordinal,Node>::CUSPSparseGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params)
    : CrsGraphBase<Ordinal,Node>(numRows,node,params)
    , isInitialized_(false)
    , isEmpty_(false)
    {
      // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
      Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
    }
    
    template <class Ordinal, class Node>
    bool CUSPSparseGraph<Ordinal,Node>::isEmpty() const
    { return isEmpty_; }
    
    template <class Ordinal, class Node>
    void CUSPSparseGraph<Ordinal,Node>::setStructure(
                        const ArrayRCP<const size_t>  &ptrs,
                        const ArrayRCP<const Ordinal> &inds)
    {
      std::string tfecfFuncName("setStructure(ptrs,inds)");
      const size_t numrows = this->getNumRows();
      const size_t numEntries = ptrs[numrows];
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          (size_t)ptrs.size() != numrows+1 
          || ptrs[0] != 0
          || (size_t)inds.size() != ptrs[numrows],
          std::runtime_error, " graph data not coherent."
      )
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          isInitialized_ == true,
          std::runtime_error, " matrix has already been initialized"
      )
      if (numrows == 0 || numEntries == 0) isEmpty_ = true;
      ptrs_ = ptrs;
      inds_ = inds;
      isInitialized_ = true;
    }
    
    template <class Ordinal, class Node>
    ArrayRCP<const size_t> CUSPSparseGraph<Ordinal,Node>::getPointers() const
    { return ptrs_; }
    
    template <class Ordinal, class Node>
    ArrayRCP<const Ordinal> CUSPSparseGraph<Ordinal,Node>::getIndices() const
    { return inds_; }
    
    template <class Ordinal, class Node>
    bool CUSPSparseGraph<Ordinal,Node>::isInitialized() const
    { return isInitialized_; }
    
    template <class Scalar, class Ordinal, class Node>
    CUSPSparseMatrix<Scalar,Ordinal,Node>::CUSPSparseMatrix(const RCP<const CUSPSparseGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params)
    : CrsMatrixBase<Scalar,Ordinal,Node>(graph,params) 
    , isInitialized_(false)
    {
      // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
      Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
    }
    
    template <class Scalar, class Ordinal, class Node>
    void CUSPSparseMatrix<Scalar,Ordinal,Node>::setValues(const ArrayRCP<const Scalar> &vals)
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
    ArrayRCP<const Scalar> CUSPSparseMatrix<Scalar,Ordinal,Node>::getValues() const
    { return vals_; }
    
    template <class Scalar, class Ordinal, class Node>
    bool CUSPSparseMatrix<Scalar,Ordinal,Node>::isInitialized() const
    {
      return isInitialized_;
    }

  };


  /** \brief Access to sparse matrix multiply and solve using the CUSP project.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class CUSPSparseOps 
  {
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
    typedef CUSPSparseOps<Scalar,Ordinal,Node> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef CUSPSparseGraph<Ordinal,Node> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef CUSPSparseMatrix<Scalar,Ordinal,Node> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c CUSPSparseOps,
    /// regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CUSPSparseOps<S2,Ordinal,Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{
  
    //! Constructor accepting and retaining a node object.
    CUSPSparseOps(const RCP<Node> &node);
  
    //! Destructor
    ~CUSPSparseOps();
  
    //@}
    //! @name Accessor routines.
    //@{ 

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Finalize a graph
    static void finalizeGraph(CUSPCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const CUSPCrsGraph<Ordinal,Node> &graph, CUSPCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(CUSPCrsGraph<Ordinal,Node> &graph, CUSPCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    void setGraphAndMatrix(const RCP<const CUSPCrsGraph<Ordinal,Node> > &graph, const RCP<const CUSPCrsMatrix<Scalar,Ordinal,Node> > &matrix);

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    CUSPSparseOps(const CUSPSparseOps& source);

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> node_;

    //! We'll fill into one and let CUSP convert into the other
    RCP<CUSP::csr_matrix<Ordinal,Scalar,CUSP::host_memory  > > hostCSR_;
    RCP<CUSP::hyb_matrix<Ordinal,Scalar,CUSP::device_memory> > devcHYB_;

    ArrayRCP<const size_t>  ptrs_;
    ArrayRCP<const Ordinal> inds_;
    ArrayRCP<const Scalar>  vals_;

    size_t numRows_, numNZ_; 
    bool isEmpty_, isInit_;

    template <class O, class S>
    static void CUSP_mult(const CUSP::hyb_matrix<Ordinal,Scalar,CUSP::device_memory> &mat, const Scalar *X, size_t ldx, Scalar *Y, size_t ldy);

    template <class O, class S> 
    static void CUSP_convert(const CUSP::crs_matrix<Ordinal,Scalar,CUSP::device_memory> *crs
                                   CUSP::hyb_matrix<Ordinal,Scalar,CUSP::device_memory> &hyb);
  };

  template<class Scalar, class Ordinal, class Node>
  CUSPSparseOps<Scalar,Ordinal,Node>::CUSPSparseOps(const RCP<Node> &node)
  : node_(node) 
  , numRows_(0)
  , numNZ_(0)
  , isEmpty_(false)
  , isInit_(false)
  { 
    // Make sure that users only specialize CUSPSparseOps for ThrustGPUNode, for now
    Teuchos::CompileTimeAssert< Teuchos::TypeTraits<Node,Kokkos::ThrustGPUNode>::is_same == false > cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node>
  CUSPSparseOps<Scalar,Ordinal,Node>::~CUSPSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> CUSPSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const {
    TEST_FOR_EXCEPTION(true, std::logic_error,
      Teuchos::typeName(*this) << "::multiply(): CUSP does not support multiple scalar types for sparse matrix-vector multiplication.");
    return;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const { 
    // beta is provided and the output multivector enjoys accumulation semantics
    TEST_FOR_EXCEPTION(true, std::logic_error,
      Teuchos::typeName(*this) << "::multiply(): CUSP does not support accumulation semantics in sparse matrix-vector multiplication.");
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    // beta is not provided and the output multivector enjoys overwrite semantics
    TEUCHOS_TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): this class does not provide support for transposed multipication. Consider manually transposing the matrix.");
    // FINISH: test for init
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node>
  template <>
  void CUSPSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                Scalar alpha, const MultiVector<Scalar,Node> &X, 
                                Scalar beta, MultiVector<Scalar,Node> &Y) const {
    // beta is provided and the output multivector enjoys accumulate semantics
    // FINISH: test for init
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): CUSP does not provide support for transposed multipication. Consider manually transposing the matrix.");
    const size_t numRHS = X.getNumCols(),
                 Xstride = X.getStride(),
                 Ystride = Y.getStride();
    ReadyBufferHelper<Node> rbh(node_);
    rbh.begin();
    const Scalar * X = rbh.template addConstBuffer<Scalar>(X.getValues());
    Scalar       * Y = rbh.template addNonConstBuffer<Scalar>(Y.getValuesNonConst());
    rbh.end();
    for (int v=0; v != numRHS; ++v) {
      CUSPSparseOps<Scalar,Ordinal,Node>::CUSP_mult<Ordinal,Scalar>(*devcHyb_,X,Y);
      X += Xstride;  
      Y += Ystride;  
    }
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */


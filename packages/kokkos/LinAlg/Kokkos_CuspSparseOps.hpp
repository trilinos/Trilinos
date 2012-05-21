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

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <cusp/hyb_matrix.h>
#include <cusp/csr_matrix.h>

namespace Kokkos {

  /** \brief Access to sparse matrix multiply and solve using the Cusp project.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class CuspSparseOps {
  public:
    //@{ 
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  ScalarType;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal OrdinalType;
    //! The Kokkos Node type.
    typedef Node    NodeType;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      // typedef whaaaa other;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      // typedef whaaaa other;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The "rebind" struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// This always specifies a specialization of \c CUSPSparseOps,
    /// regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct rebind {
      typedef CuspSparseOps<S2,Ordinal,Node> other;
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

    //! Initialize structure of matrix, using CrsGraphHostCompute
    void initializeStructure(const CrsGraphHostCompute<Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &graph);

    //! Initialize values of matrix, using CrsMatrixHostCompute
    void initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

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
    CuspSparseOps(const CuspSparseOps& source);

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> node_;

    RCP<cusp::csr_matrix<Ordinal,Scalar,cusp::host_memory  > > hostCSR_;
    RCP<cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> > devcHYB_;

    // we do this one of two ways: 
    // 1D/packed: arrays of offsets, array of ordinals, array of values.
    ArrayRCP<size_t>  begs1D_, ends1D_;
    ArrayRCP<Ordinal> inds1D_;
    ArrayRCP<Scalar>  vals1D_;
    // 2D: array of pointers
    ArrayRCP<ArrayRCP<Ordinal> > inds2D_;
    ArrayRCP<ArrayRCP<Scalar > > vals2D_;

    size_t numRows_;
    size_t numNZ_; 
    bool isEmpty_, isOpt_;
    bool indsInit_, valsInit_;

    template <class O, class S>
    static void cusp_mult(const cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> &mat, const Scalar *X, size_t ldx, Scalar *Y, size_t ldy);

    template <class O, class S> 
    static void cusp_convert(const cusp::crs_matrix<Ordinal,Scalar,cusp::device_memory> *crs
                                   cusp::hyb_matrix<Ordinal,Scalar,cusp::device_memory> &hyb);
  };

  template<class Scalar, class Ordinal, class Node>
  CuspSparseOps<Scalar,Ordinal,Node>::CuspSparseOps(const RCP<Node> &node)
  : node_(node) 
  {
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  CuspSparseOps<Scalar,Ordinal,Node>::~CuspSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> CuspSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void CuspSparseOps<Scalar,Ordinal,Node>::clear() {
    isEmpty_ = false;
    numRows_ = 0;
    numNZ_   = 0;
    indsInit_ = 0;
    valsInit_ = 0;
    begs1D_ = null;
    ends1D_ = null;
    inds1D_ = null;
    inds2D_ = null;
    //
    hostCSR_ = null;
    devcHYB_ = null;
  }

  template <class Scalar, class Ordinal, class Node>
  void CuspSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphHostCompute<Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &graph) {
  }


  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  template <class Scalar, class Ordinal, class Node>
  void CuspSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphHostCompute<Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &graph) 
  {
    TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    isOpt_ = graph.isOptimized();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.is1DStructure()) {
      isEmpty_ = false;
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t> begs, ends;
      const_cast<CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(graph).get1DStructure( inds, begs, ends );
      inds1D_ = inds;
      begs1D_ = begs;
      ends1D_ = ends;
    }
    else {
      isEmpty_  = false;
      {
        ArrayRCP<ArrayRCP<Ordinal> > inds;
        ArrayRCP<size_t>  sizes;
        const_cast<CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(graph).get2DStructure(inds,sizes);
        inds2D_     = inds;
        numEntries_ = sizes;
      }
    }
    indsInit_ = true;
  }


  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  template <class Scalar, class Ordinal, class Node>
  void CuspSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &matrix) {
    // if 2D, copy to 1D
    using Teuchos::arcp;
    TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty() || 
                       (inds2D_ != null && matrix.is1DStructure()) || (inds1D_ != null && matrix.is2DStructure()),
                       std::runtime_error, Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    ArrayRCP<const size_t>  offsets1D;
    ArrayRCP<const Ordinal> inds1D;
    ArrayRCP<const Scalar>  vals1D;
    if (!isEmpty_) 
    {
      if (matrix.is1DStructure() && isOpt_) {
        ArrayRCP<Scalar> vals;
        const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(matrix).get1DValues( vals );
        // we are awesome; nothing needs to be done
        vals1D = vals;
        begs1D = begs1D_;
        ends1D = ends1D_;
      }
      else {
        ArrayRCP<const size_t> tmpInds1D = 
        ArrayRCP<const size_t>
        if (matrix.is1DStructure()) {
          inds1D
        }
        // indices must be packed before passing to Cusp
        {
          ArrayRCP<ArrayRCP<Scalar> > vals;
          const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(matrix).get2DValues(vals);
          // convert everything (inds and vals) to 2D, temporarily
        }
      }
    }
    valsInit_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  void CuspSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &matrix) {
    TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEUCHOS_TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      ArrayRCP<Scalar> vals;
      const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,CuspSparseOps<void,Ordinal,Node> > &>(matrix).getDeviceBuffer(vals);
      pbuf_vals1D_ = vals;
    }
    valsInit_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const {
    TEST_FOR_EXCEPTION(true, std::logic_error, 
      Teuchos::typeName(*this) << "::solve(): this class does not provide support for solve.");
    return;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const { 
    TEST_FOR_EXCEPTION(true, std::logic_error,
      Teuchos::typeName(*this) << "::multiply(): Cusp does not support multiple scalar types for sparse matrix-vector multiplication.");
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void CuspSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    // beta is not provided and the output multivector enjoys overwrite semantics
    TEUCHOS_TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): this class does not provide support for transposed multipication. Consider manually transposing the matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
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
  void CuspSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                Scalar alpha, const MultiVector<Scalar,Node> &X, 
                                Scalar beta, MultiVector<Scalar,Node> &Y) const {
    // beta is provided and the output multivector enjoys accumulate semantics
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    TEST_FOR_EXCEPTION(trans != Teuchos::NO_TRANS, std::logic_error, 
      Teuchos::typeName(*this) << "::multiply(): Cusp does not provide support for transposed multipication. Consider manually transposing the matrix.");
    const size_t numRHS = X.getNumCols(),
                 Xstride = X.getStride(),
                 Ystride = Y.getStride();
    ReadyBufferHelper<Node> rbh(node_);
    rbh.begin();
    const Scalar * X = rbh.template addConstBuffer<Scalar>(X.getValues());
    Scalar       * Y = rbh.template addNonConstBuffer<Scalar>(Y.getValuesNonConst());
    rbh.end();
    for (int v=0; v != numRHS; ++v) {
      CuspSparseOps<Scalar,Ordinal,Node>::cusp_mult<Ordinal,Scalar>(*devcHyb_,X,Y);
      X += Xstride;  
      Y += Ystride;  
    }
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */


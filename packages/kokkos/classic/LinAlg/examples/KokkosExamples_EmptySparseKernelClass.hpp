/*
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
*/

#ifndef KOKKOS_EMPTY_SPARSE_KERNEL_CLASS_HPP
#define KOKKOS_EMPTY_SPARSE_KERNEL_CLASS_HPP

#include <Kokkos_CrsMatrixBase.hpp>
#include <Kokkos_CrsGraphBase.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Kokkos_DefaultSparseOps.hpp>

/// \file KokkosExamples_EmptySparseKernelClass.hpp
/// \brief A file containing a stub for a new sparse kernel provider,
///   as outlined in the \ref kokkos_crs_ops "Kokkos CRS API".

namespace KokkosExamples {

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::arcp;
  using Teuchos::ArrayView;
  using Teuchos::ParameterList;

  /// \class EmptyCrsGraph
  /// \brief An empty graph fulfilling the CrsGraph interface.
  ///
  /// This is based off Kokkos::CrsGraphBase to ease our obligations.
  template <class Node>
  class EmptyCrsGraph : public Kokkos::CrsGraphBase<int,Node> {
  public:
    EmptyCrsGraph(int numrows, int numcols, const RCP<Node> &node, const RCP<ParameterList> &params) : Kokkos::CrsGraphBase<int,Node>(numrows,numcols,node,params) {}
    ~EmptyCrsGraph() {}
    void setStructure(const ArrayRCP<const size_t>&, const ArrayRCP<const int>&) {}
  };

  /// \class EmptyCrsMatrix
  /// \brief An empty matrix fulfilling the CrsMatrix interface.
  ///
  /// This is based off Kokkos::CrsMatrixBase to ease our obligations.
  template <class Scalar, class Node>
  class EmptyCrsMatrix : public Kokkos::CrsMatrixBase<Scalar,int,Node> {
  public:
    EmptyCrsMatrix(const RCP<const EmptyCrsGraph<Node> > &graph, const RCP<ParameterList> &params) : Kokkos::CrsMatrixBase<Scalar,int,Node>(graph,params) {}
    ~EmptyCrsMatrix() {}
    void setValues(const ArrayRCP<const Scalar> &) {}
  };

  /// \class EmptySparseKernel
  /// \ingroup kokkos_crs_ops
  /// \brief Stub showing the interface that a Kokkos sparse
  ///   operations provider must implement.
  ///
  /// This class implements the Kokkos Compressed-Row Sparse API,
  /// which in turn is the interface required by the \c LocalMatOps
  /// template parameter of Tpetra::CrsMatrix.  The implementation is
  /// trivial (it does nothing), but the interface is right, so you
  /// can use it as an example for writing your own implementation of
  /// local sparse kernels.
  ///
  /// \tparam Scalar The type of entries in the sparse matrix.
  /// \tparam Node A Kokkos Node type (that implements the Kokkos Node API).
  template <class Scalar, class Node>
  class EmptySparseKernel {
  public:
    //@{
    //! @name Typedefs and structs

    /// \brief The type of entries in the sparse matrix.
    ///
    /// In your implementation, this may either be a fixed type (like
    /// <tt>double</tt>) or a template parameter of your class.
    typedef Scalar  scalar_type;
    /// \brief The type of (local) indices of the sparse matrix.
    ///
    /// This is \c void only because this is a stub implementation.
    /// In a real implementation, ordinal_type would normally either be
    /// a fixed type (like \c int) or a template parameter of your
    /// class.
    typedef void ordinal_type;
    //! The Kokos Node type.
    typedef Node    node_type;
    //! The type of this object: <tt>typeof(*this)</tt>
    typedef EmptySparseKernel<Scalar,Node> sparse_ops_type;

    //! Typedef for local graph class
    template <class O, class N>
    struct graph {
      typedef EmptyCrsGraph<N> graph_type;
    };

    //! Typedef for local matrix class
    template <class S, class O, class N>
    struct matrix {
      typedef EmptyCrsMatrix<S,N> matrix_type;
    };

    /// \brief Rebind struct, for specifying type information for a different scalar.
    ///
    /// This typedef lets you tell us where to find sparse kernels for
    /// sparse matrices with entries of scalar type T.  T may be
    /// different than scalar_type.
    ///
    /// The point of this typedef is to provide a "fall-back"
    /// implementation as a function of the scalar type T.  This is
    /// because you might have specific kernels for some scalar types
    /// T (such as <tt>float</tt> or <tt>double</tt>), but not for
    /// other types (such as extended-precision types).  For example,
    /// you might want to use a third-party library's kernels, but the
    /// library only implements the required kernels for scalar types
    /// <tt>float</tt> and <tt>double</tt>.  Thus, you can provide a
    /// "fall-back" implementation for scalar types not supported by
    /// the third-party library.  We do this for Kokkos::MklSparseOps,
    /// for example.  Alternately, some scalar types T (especially
    /// those requiring dynamic memory allocation) might not work
    /// correctly or efficiently on certain Kokkos Node types
    /// (especially GPU Node types).  In that case, you could use the
    /// fall-back implementation to generate either a compile-time or
    /// a run-time error for unsupported scalar types.
    template <class T>
    struct bind_scalar {
      typedef EmptySparseKernel<T,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a Kokkos Node instance.
    EmptySparseKernel(const RCP<Node> & node) : node_(node) {}

    //! Destructor.
    ~EmptySparseKernel() {}

    //@}
    //! @name Accessor routines.
    //@{

    /// \brief Kokkos Node accessor.
    ///
    /// Return the Kokkos Node instance of type this::node_type given
    /// to the constructor.
    RCP<Node> getNode() const {return node_;}

    //@}
    //! @name Initialization of structure
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
    {
      return Kokkos::details::DefaultCRSAllocator::template allocRowPtrs<size_t,Node>(node,numEntriesPerRow);
    }

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
    {
      return Kokkos::details::DefaultCRSAllocator::template allocStorage<T,int,Node>(node,rowPtrs);
    }

    //! Finalize a graph
    inline static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, EmptyCrsGraph<Node> &graph, const RCP<ParameterList> &params) {}

    //! Finalize the matrix of an already-finalized graph.
    inline static void finalizeMatrix(const EmptyCrsGraph<Node> &graph, EmptyCrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params) {}

    //! Finalize a graph and a matrix.
    inline static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, EmptyCrsGraph<Node> &graph, EmptyCrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params) {}

    /** \brief Initialize the kernels with the graph and matrix.

        This is the mechanism by which the user specifies the
        structure and values for the sparse matrix ops.

        setGraphAndMatrix() must be called before calling multiply() or solve().

        After initializeStructure() completes, the caller is
        responsible for deciding what to do with the graph and matrix objects.
        Since your implementation may choose just
        to view the original CrsGraph data instead of making a deep
        copy, callers should not change the Kokkos::CrsGraph after
        calling this method.
      */
    void setGraphAndMatrix(const RCP<const EmptyCrsGraph<Node> >         &graph,
                           const RCP<const EmptyCrsMatrix<Scalar,Node> > &node) {};

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
              const Kokkos::MultiVector<DomainScalar,Node> &X,
              Kokkos::MultiVector<RangeScalar,Node> &Y) const
    {}

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
              const Kokkos::MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              Kokkos::MultiVector<RangeScalar,Node> &Y) const
    {}

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
           const Kokkos::MultiVector<DomainScalar,Node> &Y,
           Kokkos::MultiVector<RangeScalar,Node> &X) const
    {}

    /// \brief Gauss-Seidel or SOR on \f$B = A X\f$.
    ///
    /// Apply a forward or backward sweep of Gauss-Seidel or
    /// Successive Over-Relaxation (SOR) to the linear system(s) \f$B
    /// = A X\f$.  For Gauss-Seidel, set the damping factor \c omega
    /// to 1.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in B.
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector B.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param omega [in] SOR damping factor.  omega = 1 results in
    ///   Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward or Backward.
    ///   If you want a symmetric sweep, call this method twice, first
    ///   with direction = Forward then with direction = Backward.
    ///
    /// \note We don't include a separate "Symmetric" direction mode
    ///   in order to avoid confusion when using this method to
    ///   implement "hybrid" Jacobi + symmetric (Gauss-Seidel or SOR)
    ///   for a matrix distributed over multiple processes.  ("Hybrid"
    ///   means "Gauss-Seidel or SOR within the process, Jacobi
    ///   outside.")  In that case, interprocess communication (a
    ///   boundary exchange) must occur before both the forward sweep
    ///   and the backward sweep, so we would need to invoke the
    ///   kernel once per sweep direction anyway.
    template <class DomainScalar, class RangeScalar>
    void
    gaussSeidel (const Kokkos::MultiVector<DomainScalar,Node> &B,
                 Kokkos::MultiVector<RangeScalar,Node> &X,
                 const Kokkos::MultiVector<Scalar,Node> &D,
                 const RangeScalar& dampingFactor,
                 const Kokkos::ESweepDirection direction) const
    {}

    //@}
  protected:
    //! Copy constructor (protected and unimplemented)
    EmptySparseKernel(const EmptySparseKernel& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;
  };

  /** \example EmptySparseKernelDriver.cpp
   *
   * This example demonstrates the basic use case for a sparse kernel
   * provider. It also verifies that the stub builds correctly.
   */

} // end of namespace KokkosExamples

#endif

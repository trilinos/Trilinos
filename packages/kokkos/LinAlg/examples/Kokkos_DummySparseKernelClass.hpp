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

#ifndef KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP
#define KOKKOS_DUMMY_SPARSE_KERNEL_CLASS_HPP

#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_CrsGraph.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_BLAS_types.hpp>

/// \file Kokkos_DummySparseKernelClass.hpp
/// \brief A file containing a stub for a new sparse kernel provider,
///   as outlined in the \ref kokkos_crs_ops "Kokkos CRS API".

namespace KokkosExamples {

  /// \class DummySparseKernel
  /// \ingroup kokkos_crs_ops
  /// \brief Stub showing the interface that a Kokkos sparse
  ///   operations provider must implement.
  ///
  /// This class implements the Kokkos Compressed-Row Sparse API,
  /// which in turn is the interface required by the \c LocalMatOps
  /// template parameter of \c Tpetra::CrsMatrix.  The implementation
  /// is trivial (it does nothing), but the interface is right, so you
  /// can use it as an example for writing your own implementation of
  /// local sparse kernels.
  ///
  /// \tparam Node A Kokkos Node type (that implements the Kokkos Node API).
  template <class Node>
  class DummySparseKernel {
  public:
    //@{
    //! @name Typedefs and structs

    /// \brief The type of entries of the sparse matrix.
    ///
    /// This is \c void only because this is a stub implementation.
    /// In a real implementation, ScalarType would normally either be
    /// a fixed type (like \c double) or a template parameter of your
    /// class.
    typedef void  ScalarType;
    /// \brief The type of (local) indices of the sparse matrix.
    ///
    /// This is \c void only because this is a stub implementation.
    /// In a real implementation, OrdinalType would normally either be
    /// a fixed type (like \c int) or a template parameter of your
    /// class.
    typedef void OrdinalType;
    //! The Kokos Node type.
    typedef Node    NodeType;
    //! The type of this object: \c typeof(*this)
    typedef DummySparseKernel<Node> ThisType;

    /// \brief Rebind struct, for specifying type information for a different scalar.
    ///
    /// This typedef lets you tell us where to find sparse kernels for
    /// sparse matrices with entries of scalar type T.  T may be
    /// different than ScalarType.
    ///
    /// One point of this typedef is that sometimes you may have
    /// noptimized kernels for some scalar types T (such as float or
    /// double), but not for other types (such as extended-precision
    /// types).  Some scalar types T (especially those requiring
    /// dynamic memory allocation) might not work correctly or
    /// efficiently on certain Kokkos Node types (especially GPU Node
    /// types).  This typedef lets you provide a "fall-back"
    /// implementation of sparse kernels.
    template <class T>
    struct rebind {
      typedef DummySparseKernel<Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a Kokkos Node instance.
    DummySparseKernel(const Teuchos::RCP<Node> & node) : node_(node) {}

    //! Destructor.
    ~DummySparseKernel() {}

    //@}
    //! @name Accessor routines.
    //@{

    /// \brief Kokkos Node accessor.
    ///
    /// Return the Kokkos Node instance of type this::NodeType given
    /// to the constructor.
    Teuchos::RCP<Node> getNode() const {return node_;}

    //@}
    //! @name Initialization of structure
    //@{

    /** \brief Initialize the structure of the sparse matrix.

        This is the mechanism by which the user specifies the
        structure for the sparse matrix.  The structure always come
        via a Kokkos::CrsGraph<O,N,SO> instance, where

        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        You must first provide the graph structure via
        initializeStructure() before you can provide the matrix values
        to initializeValues(). After calling initializeStructure(),
        the clear() method must be called before you can call
        initializeStructure() again.

        In general, both initializeStructure() and initializeValues()
        must be called before calling multiply() or solve().

        Note that your implementation of this method is not obligated
        to represent the sparse matrix's structure in the same way
        that Kokkos::CrsGraph does.  That is just an entry format.
        Your implementation may choose either to make a deep copy of
        the input data (and possibly change the storage format), or
        simply to view it (by copying the ArrayRCPs).

        After initializeStructure() completes, the caller is
        responsible for deciding what to do with the original
        Kokkos::CrsGraph.  Since your implementation may choose just
        to view the original CrsGraph data instead of making a deep
        copy, callers should not change the Kokkos::CrsGraph after
        calling this method, unless they first call clear().
      */
    template<class Ordinal>
    void
    initializeStructure (const Kokkos::CrsGraph<Ordinal,Node,DummySparseKernel<Node> > & /* graph */) {}

    /** \brief Initialize the values of the sparse matrix.

        This is the mechanism by which the user specifies the values
        for the sparse matrix.  The values always come via a
        Kokkos::CrsMatrix<S,O,N,SO> instance, where

        - \c S is the scalar type this::ScalarType
        - \c O is the ordinal type this::OrdinalType,
        - \c N is the node type this::NodeType, and
        - \c SO is the sparse op type this::ThisType.

        You must first provide the graph structure via
        initializeStructure() before you can provide the matrix values
        to initializeValues().  You may provide the matrix values
        repeatedly via multiple calls to initializeValues() without
        needing to call clear() in between.

        In general, both initializeStructure() and initializeValues()
        must be called before calling multiply() or solve().

        Note that your implementation of this method is not obligated
        to represent the sparse matrix's entries in the same way that
        Kokkos::CrsMatrix does.  That is just an entry format.  Your
        implementation may choose either to make a deep copy of the
        input data (and possibly change the storage format), or simply
        to view it (by copying the ArrayRCPs).  For example, you might
        choose to copy the graph structure and values from their input
        compressed sparse row format into jagged diagonal storage.

        After initializeValues() completes, the caller is responsible
        for deciding what to do with the original Kokkos::CrsMatrix.
        Since your implementation may choose just to view the original
        CrsGraph data instead of making a deep copy, callers should
        not change the Kokkos::CrsMatrix after calling this method,
        unless they first call clear().
    */
    template<class Scalar, class Ordinal>
    void
    initializeValues (const Kokkos::CrsMatrix<Scalar,Ordinal,Node,DummySparseKernel<Node> > & /* matrix */) {}

    /// \brief Clear the graph and matrix data.
    ///
    /// If you have already called initializeStructure(), you must
    /// first call clear() before you may call initializeStructure()
    /// again.
    ///
    /// After calling clear(), no significant data (including
    /// persisting references) will be preserved, save for the pointer
    /// to the Kokkos Node instance.  You must then call
    /// initializeStructure() and initializeValues() (in that order)
    /// again before you may call multiply() or solve().
    void clear () {}

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
           Teuchos::EUplo uplo,
           Teuchos::EDiag diag,
           const Kokkos::MultiVector<DomainScalar,Node> &Y,
           Kokkos::MultiVector<RangeScalar,Node> &X) const
    {}

    //@}
  protected:
    //! Copy constructor (protected and unimplemented)
    DummySparseKernel(const DummySparseKernel& source);

    //! The Kokkos Node instance given to this object's constructor.
    Teuchos::RCP<Node> node_;
  };

  /** \example DummySparseKernelDriver.cpp
   *
   * This example demonstrates the basic use case for a sparse kernel
   * provider. It also verifies that the stub builds correctly.
   */

} // end of namespace KokkosExamples

#endif

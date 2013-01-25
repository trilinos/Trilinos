// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include <Teuchos_TimeMonitor.hpp>


/*! \file Tpetra_CrsMatrixMultiplyOp_decl.hpp

    The declarations for the class Tpetra::CrsMatrixMultiplyOp and related non-member constructors.
 */

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration
  template <class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class CrsMatrix;
#endif

  /// \brief A class for wrapping a CrsMatrix multiply in a Operator.
  ///
  /// \note Most Tpetra users do not need to use this class.  It will
  ///   be useful to Tpetra users who want to do mixed-precision
  ///   sparse matrix-vector multiply, where the sparse matrix's
  ///   entries have a different precision than that of the input and
  ///   output vectors.  If your sparse matrix and vectors have the
  ///   same type of entries, then you don't need to use this class.
  ///
  /// This class makes a <tt>CrsMatrix<MatScalar, ...></tt> "look
  /// like" an <tt>Operator<Scalar, ...></tt>, where
  /// <tt>MatScalar</tt> and <tt>Scalar</tt> may be different types.
  /// It does so by working around a limitation of C++, namely that
  /// template methods of a class can't be virtual.
  ///
  /// Here is a detailed description of how the language issue relates
  /// to CrsMatrix.  If you call the <tt>apply</tt> method of
  /// CrsMatrix, you will always get the version that takes a
  /// <tt>MultiVector<Scalar, ...></tt> input and produces a
  /// <tt>MultiVector<Scalar, ...></tt> output.  CrsMatrix actually
  /// implements a a templated sparse matrix-vector multiply operation
  /// (its <tt>localMultiply</tt> method).  It is templated on the
  /// scalar types of its input and output multivectors
  /// (<tt>DomainScalar</tt> resp. <tt>RangeScalar</tt>).  However,
  /// Operator can't access this templated mat-vec method.  This is
  /// because Operator::apply is virtual, and therefore cannot have a
  /// template parameter for the <tt>Scalar</tt> type of the
  /// MultiVector input and output.
  ///
  /// Users who want to access the templated sparse mat-vec in
  /// CrsMatrix through the Operator interface may wrap the CrsMatrix
  /// in an instance of this class.  This class implements an Operator
  /// that takes <tt>MultiVector<Scalar, ...></tt> input and output,
  /// but the CrsMatrix may contain any desired type
  /// <tt>MatScalar</tt>.  The type <tt>MatScalar</tt> may differ from
  /// the <tt>Scalar</tt> type of the MultiVector input and output.
  /// That works around the "no virtual template methods" issue for
  /// input and output multivectors of the same type.
  ///
  /// Interestingly enough, CrsMatrix implements its <tt>apply</tt>
  /// method using an instance of this class with <tt>Scalar ==
  /// MatScalar</tt>.  CrsMatrix does not actually contain an
  /// implementation of "nonlocal" (distributed over multiple MPI
  /// processes) mat-vec; its <tt>apply</tt> defers the nonlocal part
  /// to this class' apply() method.  The same is true for the
  /// gaussSeidel() method.
  ///
  /// \tparam Scalar The type of the entries of the input and output
  ///   MultiVector of the apply() method.  Same as the first template
  ///   parameter of Operator.
  ///
  /// \tparam MatScalar The type of the entries of the CrsMatrix; the
  ///   first template parameter of CrsMatrix.
  ///
  /// \tparam LocalOrdinal The second template parameter of CrsMatrix
  ///   and Operator.
  ///
  /// \tparam GlobalOrdinal The third template parameter of CrsMatrix
  ///   and Operator.
  ///
  /// \tparam Node The fourth template parameter of CrsMatrix and
  ///   Operator.
  ///
  /// \tparam LocalMatOps The fifth template parameter of CrsMatrix.
  ///   (Operator only takes four template parameters.)
  template <class Scalar,
            class MatScalar = Scalar,
            class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename Kokkos::DefaultKernels<MatScalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrixMultiplyOp :
    public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  {
  public:
    //! @name Constructor and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param A [in] The CrsMatrix to wrap as an
    ///   <tt>Operator<Scalar, ...></tt>.
    CrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

    //! Destructor
    virtual ~CrsMatrixMultiplyOp();

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute <tt>Y = beta*Y + alpha*Op(A)*X</tt>, where
    ///   <tt>Op(A)</tt> is either A, \f$A^T\f$, or \f$A^H\f$.
    ///
    /// This method calls the underlying CrsMatrix object's
    /// localMultiply<Scalar,Scalar>() method.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    /// \brief "Hybrid" Jacobi + (Gauss-Seidel or SOR) on \f$B = A X\f$.
    ///
    /// "Hybrid" means Jacobi for interprocess communication, but
    /// Successive Over-Relaxation (SOR) or Gauss-Seidel for
    /// intraprocess computation.  Gauss-Seidel is a special case of
    /// SOR, where the damping factor is one.
    ///
    /// The Forward or Backward sweep directions have their usual SOR
    /// meaning within the process.  Interprocess communication occurs
    /// once before the sweep, as it would in Jacobi.
    ///
    /// The Symmetric sweep direction means first Forward, then
    /// Backward.  Before each sweep is an interprocess communication,
    /// as in Jacobi.  Thus, Symmetric results in two interprocess
    /// communication steps.
    ///
    /// \param B [in] Right-hand side(s), in the range Map of the
    ///   matrix.
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).  This must be a domain Map view of
    ///   a column Map multivector.
    /// \param D [in] Inverse of diagonal entries of the matrix A,
    ///   in the row Map of the matrix.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its Backward
    ///   sweep) as one.
    ///
    /// \pre Domain, range, and row Maps of the sparse matrix are all
    ///   the same.  (The domain and range Maps must be the same
    ///   because this kernel overwrites its input.  The row Map must
    ///   be the same because the kernel uses the same local indices
    ///   for the rows of the sparse matrix, and for the rows of the
    ///   input / output multivector.)
    ///
    /// \pre No other argument aliases X.
    void
    gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                 const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                 const Scalar& dampingFactor,
                 const ESweepDirection direction,
                 const int numSweeps) const;

    /// \brief Version of gaussSeidel(), with fewer requirements on X.
    ///
    /// This method is just like gaussSeidel(), except that X need
    /// only be in the domain Map.  This method does not require that
    /// X be a domain Map view of a column Map multivector.  As a
    /// result, this method must copy X into a domain Map multivector
    /// before operating on it.
    ///
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param B [in] Right-hand side(s), in the range Map.
    /// \param D [in] Inverse of diagonal entries of the matrix,
    ///   in the row Map.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its
    ///   Backward sweep) as one.
    ///
    /// \pre Domain, range, and row Maps of the sparse matrix are
    ///   all the same.
    /// \pre No other argument aliases X.
    void
    gaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                     const Scalar& dampingFactor,
                     const ESweepDirection direction,
                     const int numSweeps) const;

    /// \brief Whether this Operator's apply() method can apply the
    ///   transpose or conjugate transpose.
    ///
    /// This is always true, since it is true for the CrsMatrix that
    /// this object wraps.
    bool hasTransposeApply() const;

    //! The domain Map of this Operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

    //! The range Map of this Operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

    //@}

  protected:
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix object.
    const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > matrix_;

    /// \brief Column Map MultiVector used in apply().
    ///
    /// This is a column Map MultiVector.  It is used as the target of
    /// the forward mode Import operation (if necessary) in
    /// applyNonTranspose(), and the source of the reverse mode Export
    /// operation (if necessary) in applyTranspose().  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the apply() method
    /// is const, yet the method needs to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importMV_;

    /// \brief Row Map MultiVector used in apply().
    ///
    /// This is a row Map MultiVector.  It is uses as the source of
    /// the forward mode Export operation (if necessary) in
    /// applyNonTranspose(), and the target of the reverse mode Import
    /// operation (if necessary) in applyTranspose().  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the apply() method
    /// is const, yet the method needs to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > exportMV_;

#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    Teuchos::RCP<Teuchos::Time> importTimer_, exportTimer_;
#endif

    /// \brief Apply the transpose or conjugate transpose of the
    ///   matrix to X, producing Y.
    void
    applyTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                    Scalar alpha,
                    Scalar beta) const;

    //! Apply the matrix (not its transpose) to X, producing Y.
    void
    applyNonTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                       MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                       Scalar alpha,
                       Scalar beta) const;

  private:
    /// \brief Create a (or fetch a cached) column Map MultiVector.
    ///
    /// \param X_domainMap [in] A domain Map Multivector.  The
    ///   returned MultiVector, if nonnull, will have the same number
    ///   of columns as Y_domainMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the domain Map and the
    /// column Map are the same (so that normally we wouldn't need the
    /// column Map MultiVector), but the following (for example)
    /// holds:
    ///
    /// 1. The kernel needs a constant stride input MultiVector, but
    ///    the given input MultiVector is not constant stride.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getColumnMapMultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X_domainMap,
                             const bool force = false) const;

    /// \brief Create a (or fetch a cached) row Map MultiVector.
    ///
    /// \param Y_rangeMap [in] A range Map Multivector.  The returned
    ///   MultiVector, if nonnull, will have the same number of
    ///   columns as Y_rangeMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the range Map and the
    /// row Map are the same (so that normally we wouldn't need the
    /// row Map MultiVector), but one of the following holds:
    ///
    /// 1. The kernel needs a constant stride output MultiVector,
    ///    but the given output MultiVector is not constant stride.
    ///
    /// 2. The kernel does not permit aliasing of its input and output
    ///    MultiVector arguments, but they do alias each other.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getRowMapMultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y_rangeMap,
                          const bool force = false) const;
  };

  /// \brief Non-member function to create a CrsMatrixMultiplyOp.
  /// \relatesalso CrsMatrixMultiplyOp
  ///
  /// The function has the same template parameters of CrsMatrixMultiplyOp.
  ///
  /// \param A [in] The CrsMatrix instance to wrap in an CrsMatrixMultiplyOp.
  /// \return The CrsMatrixMultiplyOp wrapper for the given CrsMatrix.
  template <class Scalar,
            class MatScalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  Teuchos::RCP< CrsMatrixMultiplyOp<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
  createCrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP

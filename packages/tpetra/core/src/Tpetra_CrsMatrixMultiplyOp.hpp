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

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_HPP

/// \file Tpetra_CrsMatrixMultiplyOp.hpp
///
/// Declaration and definition of Tpetra::CrsMatrixMultiplyOp and its
/// nonmember constructor Tpetra::createCrsMatrixMultiplyOp.

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"

namespace Tpetra {

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
  /// \tparam LocalOrdinal The type of the local indices of the
  ///   CrsMatrix; the second template parameter of CrsMatrix and
  ///   Operator.
  ///
  /// \tparam GlobalOrdinal The type of the global indices of the
  ///   CrsMatrix; the third template parameter of CrsMatrix and
  ///   Operator.
  ///
  /// \tparam Node The fourth template parameter of CrsMatrix and
  ///   Operator.
  template <class Scalar,
            class MatScalar = Scalar,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
  class CrsMatrixMultiplyOp :
    public Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    //! The specialization of CrsMatrix which this class wraps.
    typedef CrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    //! The specialization of Map which this class uses.
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    //! @name Constructor and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param A [in] The CrsMatrix to wrap as an
    ///   <tt>Operator<Scalar, ...></tt>.
    CrsMatrixMultiplyOp (const Teuchos::RCP<const crs_matrix_type>& A) :
      matrix_ (A)
    {
      // we don't require that A is fill complete; we will query for the
      // importer/exporter at apply()-time
    }

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~CrsMatrixMultiplyOp () {}

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
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! matrix_->isFillComplete (), std::runtime_error,
         Teuchos::typeName (*this) << "::apply(): underlying matrix is not fill-complete.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
         Teuchos::typeName (*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
         Teuchos::typeName (*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
      if (mode == Teuchos::NO_TRANS) {
        applyNonTranspose (X, Y, alpha, beta);
      }
      else {
        applyTranspose (X, Y, mode, alpha, beta);
      }
    }

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
                 const int numSweeps) const
    {
      using Teuchos::null;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      using Teuchos::rcp_const_cast;
      typedef Scalar OS;
      typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
      typedef Export<LocalOrdinal, GlobalOrdinal, Node> export_type;
      typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
      typedef MultiVector<OS, LocalOrdinal, GlobalOrdinal, Node> OSMV;

      TEUCHOS_TEST_FOR_EXCEPTION
        (numSweeps < 0, std::invalid_argument,
         "gaussSeidel: The number of sweeps must be nonnegative, "
         "but you provided numSweeps = " << numSweeps << " < 0.");

      // Translate from global to local sweep direction.
      // While doing this, validate the input.
      ESweepDirection localDirection;
      if (direction == Forward) {
        localDirection = Forward;
      }
      else if (direction == Backward) {
        localDirection = Backward;
      }
      else if (direction == Symmetric) {
        // We'll control local sweep direction manually.
        localDirection = Forward;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument,
           "gaussSeidel: The 'direction' enum does not have any of its valid "
           "values: Forward, Backward, or Symmetric.");
      }

      if (numSweeps == 0) {
        return; // Nothing to do.
      }

      // We don't need the Export object because this method assumes
      // that the row, domain, and range Maps are the same.  We do need
      // the Import object, if there is one, though.
      RCP<const import_type> importer = matrix_->getGraph()->getImporter();
      RCP<const export_type> exporter = matrix_->getGraph()->getExporter();
      TEUCHOS_TEST_FOR_EXCEPTION
        (! exporter.is_null (), std::runtime_error,
         "Tpetra's gaussSeidel implementation requires that the row, domain, "
         "and range Maps be the same.  This cannot be the case, because the "
         "matrix has a nontrivial Export object.");

      RCP<const map_type> domainMap = matrix_->getDomainMap ();
      RCP<const map_type> rangeMap = matrix_->getRangeMap ();
      RCP<const map_type> rowMap = matrix_->getGraph ()->getRowMap ();
      RCP<const map_type> colMap = matrix_->getGraph ()->getColMap ();

      const bool debug = ::Tpetra::Details::Behavior::debug ();
      if (debug) {
        // The relation 'isSameAs' is transitive.  It's also a
        // collective, so we don't have to do a "shared" test for
        // exception (i.e., a global reduction on the test value).
        TEUCHOS_TEST_FOR_EXCEPTION
          (! X.getMap ()->isSameAs (*domainMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidel requires that the input "
           "multivector X be in the domain Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! B.getMap ()->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidel requires that the input "
           "B be in the range Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! D.getMap ()->isSameAs (*rowMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidel requires that the input "
           "D be in the row Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! rowMap->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidel requires that the row Map and the "
           "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! domainMap->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidel requires that the domain Map and "
           "the range Map of the matrix be the same.");
      }

      // If B is not constant stride, copy it into a constant stride
      // multivector.  We'l handle the right-hand side B first and deal
      // with X right before the sweeps, to improve locality of the
      // first sweep.  (If the problem is small enough, then that will
      // hopefully keep more of the entries of X in cache.  This
      // optimizes for the typical case of a small number of sweeps.)
      RCP<const OSMV> B_in;
      if (B.isConstantStride()) {
        B_in = rcpFromRef (B);
      }
      else {
        // The range Map and row Map are the same in this case, so we
        // can use the (possibly cached) row Map multivector to store a
        // constant stride copy of B.  We don't have to copy back, since
        // Gauss-Seidel won't modify B.
        RCP<OSMV> B_in_nonconst = getRowMapMultiVector (B, true);
        deep_copy (*B_in_nonconst, B);
        B_in = rcp_const_cast<const OSMV> (B_in_nonconst);

        TPETRA_EFFICIENCY_WARNING
          (! B.isConstantStride (), std::runtime_error,
           "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
           "requires that X and B both have constant stride.  Since B does not "
           "have constant stride, we had to make a copy.  This is a limitation of "
           "the current implementation and not your fault, but we still report it "
           "as an efficiency warning for your information.");
      }

      // If X is not constant stride, copy it into a constant stride
      // multivector.  Also, make the column Map multivector X_colMap,
      // and its domain Map view X_domainMap.  (X actually must be a
      // domain Map view of a column Map multivector; exploit this, if X
      // has constant stride.)

      RCP<OSMV> X_domainMap;
      RCP<OSMV> X_colMap;
      bool copiedInput = false;

      if (importer.is_null ()) { // Domain and column Maps are the same.
        if (X.isConstantStride ()) {
          X_domainMap = rcpFromRef (X);
          X_colMap = X_domainMap;
          copiedInput = false;
        }
        else {
          // Get a temporary column Map multivector, make a domain Map
          // view of it, and copy X into the domain Map view.  We have
          // to copy here because we won't be doing Import operations.
          X_colMap = getColumnMapMultiVector (X, true);
          X_domainMap = X_colMap; // Domain and column Maps are the same.
          deep_copy (*X_domainMap, X); // Copy X into the domain Map view.
          copiedInput = true;
          TPETRA_EFFICIENCY_WARNING
            (! X.isConstantStride (), std::runtime_error,
             "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
             "requires that X and B both have constant stride.  Since X does not "
             "have constant stride, we had to make a copy.  This is a limitation of "
             "the current implementation and not your fault, but we still report it "
             "as an efficiency warning for your information.");
        }
      }
      else { // We will be doing Import operations in the sweeps.
        if (X.isConstantStride ()) {
          X_domainMap = rcpFromRef (X);
          // This kernel assumes that X is a domain Map view of a
          // column Map multivector.  We will only check if this is
          // valid in debug mode.
          X_colMap = X_domainMap->offsetViewNonConst (colMap, 0);

          // Do the first Import for the first sweep.  This simplifies
          // the logic in the sweeps.
          X_colMap->doImport (X, *importer, INSERT);
          copiedInput = false;
        }
        else {
          // Get a temporary column Map multivector X_colMap, and make a
          // domain Map view X_domainMap of it.  Instead of copying, we
          // do an Import from X into X_domainMap.  This saves us a
          // copy, since the Import has to copy the data anyway.
          X_colMap = getColumnMapMultiVector (X, true);
          X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);
          X_colMap->doImport (X, *importer, INSERT);
          copiedInput = true;
          TPETRA_EFFICIENCY_WARNING
            (! X.isConstantStride (), std::runtime_error,
             "gaussSeidel: The current implementation of the Gauss-Seidel kernel "
             "requires that X and B both have constant stride.  Since X does not "
             "have constant stride, we had to make a copy.  This is a limitation of "
             "the current implementation and not your fault, but we still report it "
             "as an efficiency warning for your information.");
        }
      }

      for (int sweep = 0; sweep < numSweeps; ++sweep) {
        if (! importer.is_null () && sweep > 0) {
          // We already did the first Import for the zeroth sweep.
          X_colMap->doImport (*X_domainMap, *importer, INSERT);
        }

        // Do local Gauss-Seidel.
        if (direction != Symmetric) {
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     localDirection);
        }
        else { // direction == Symmetric
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     Forward);
          // Communicate again before the Backward sweep.
          if (! importer.is_null ()) {
            X_colMap->doImport (*X_domainMap, *importer, INSERT);
          }
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     Backward);
        }
      }

      if (copiedInput) {
        deep_copy (X, *X_domainMap); // Copy back: X_domainMap -> X.
      }
    }

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
                     const int numSweeps) const
    {
      using Teuchos::null;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      using Teuchos::rcp_const_cast;
      typedef Scalar OS;
      typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
      typedef Export<LocalOrdinal, GlobalOrdinal, Node> export_type;
      typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
      typedef MultiVector<OS, LocalOrdinal, GlobalOrdinal, Node> OSMV;

      TEUCHOS_TEST_FOR_EXCEPTION
        (numSweeps < 0, std::invalid_argument,
         "gaussSeidelCopy: The number of sweeps must be nonnegative, "
         "but you provided numSweeps = " << numSweeps << " < 0.");

      // Translate from global to local sweep direction.
      // While doing this, validate the input.
      ESweepDirection localDirection;
      if (direction == Forward) {
        localDirection = Forward;
      }
      else if (direction == Backward) {
        localDirection = Backward;
      }
      else if (direction == Symmetric) {
        // We'll control local sweep direction manually.
        localDirection = Forward;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument,
           "gaussSeidelCopy: The 'direction' enum does not have any of its "
           "valid values: Forward, Backward, or Symmetric.");
      }

      if (numSweeps == 0) {
        return;
      }

      RCP<const import_type> importer = matrix_->getGraph()->getImporter();
      RCP<const export_type> exporter = matrix_->getGraph()->getExporter();
      TEUCHOS_TEST_FOR_EXCEPTION
        (! exporter.is_null (),
         std::runtime_error,
         "Tpetra's gaussSeidelCopy implementation requires that the row, domain, "
         "and range Maps be the same.  This cannot be the case, because the "
         "matrix has a nontrivial Export object.");

      RCP<const map_type> domainMap = matrix_->getDomainMap ();
      RCP<const map_type> rangeMap = matrix_->getRangeMap ();
      RCP<const map_type> rowMap = matrix_->getGraph ()->getRowMap ();
      RCP<const map_type> colMap = matrix_->getGraph ()->getColMap ();

      const bool debug = ::Tpetra::Details::Behavior::debug ();
      if (debug) {
        // The relation 'isSameAs' is transitive.  It's also a
        // collective, so we don't have to do a "shared" test for
        // exception (i.e., a global reduction on the test value).
        TEUCHOS_TEST_FOR_EXCEPTION
          (! X.getMap ()->isSameAs (*domainMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
           "multivector X be in the domain Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! B.getMap ()->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
           "B be in the range Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! D.getMap ()->isSameAs (*rowMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidelCopy requires that the input "
           "D be in the row Map of the matrix.");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! rowMap->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidelCopy requires that the row Map and the "
           "range Map be the same (in the sense of Tpetra::Map::isSameAs).");
        TEUCHOS_TEST_FOR_EXCEPTION
          (! domainMap->isSameAs (*rangeMap), std::runtime_error,
           "Tpetra::CrsMatrix::gaussSeidelCopy requires that the domain Map and "
           "the range Map of the matrix be the same.");
      }

      // Fetch a (possibly cached) temporary column Map multivector
      // X_colMap, and a domain Map view X_domainMap of it.  Both have
      // constant stride by construction.  We know that the domain Map
      // must include the column Map, because our Gauss-Seidel kernel
      // requires that the row Map, domain Map, and range Map are all
      // the same, and that each process owns all of its own diagonal
      // entries of the matrix.

      RCP<OSMV> X_colMap;
      RCP<OSMV> X_domainMap;
      bool copyBackOutput = false;
      if (importer.is_null ()) {
        if (X.isConstantStride ()) {
          X_colMap = rcpFromRef (X);
          X_domainMap = rcpFromRef (X);
          // No need to copy back to X at end.
        }
        else { // We must copy X into a constant stride multivector.
          // Just use the cached column Map multivector for that.
          X_colMap = getColumnMapMultiVector (X, true);
          // X_domainMap is always a domain Map view of the column Map
          // multivector.  In this case, the domain and column Maps are
          // the same, so X_domainMap _is_ X_colMap.
          X_domainMap = X_colMap;
          deep_copy (*X_domainMap, X); // Copy X into constant stride multivector
          copyBackOutput = true; // Don't forget to copy back at end.
          TPETRA_EFFICIENCY_WARNING
            (! X.isConstantStride (), std::runtime_error,
             "gaussSeidelCopy: The current implementation of the Gauss-Seidel "
             "kernel requires that X and B both have constant stride.  Since X "
             "does not have constant stride, we had to make a copy.  This is a "
             "limitation of the current implementation and not your fault, but we "
             "still report it as an efficiency warning for your information.");
        }
      }
      else { // Column Map and domain Map are _not_ the same.
        X_colMap = getColumnMapMultiVector (X);
        X_domainMap = X_colMap->offsetViewNonConst (domainMap, 0);

        // We could just copy X into X_domainMap.  However, that wastes
        // a copy, because the Import also does a copy (plus
        // communication).  Since the typical use case for Gauss-Seidel
        // is a small number of sweeps (2 is typical), we don't want to
        // waste that copy.  Thus, we do the Import here, and skip the
        // first Import in the first sweep.  Importing directly from X
        // effects the copy into X_domainMap (which is a view of
        // X_colMap).
        X_colMap->doImport (X, *importer, INSERT);

        copyBackOutput = true; // Don't forget to copy back at end.
      }

      // The Gauss-Seidel / SOR kernel expects multivectors of constant
      // stride.  X_colMap is by construction, but B might not be.  If
      // it's not, we have to make a copy.
      RCP<const OSMV> B_in;
      if (B.isConstantStride ()) {
        B_in = rcpFromRef (B);
      }
      else {
        // Range Map and row Map are the same in this case, so we can
        // use the cached row Map multivector to store a constant stride
        // copy of B.
        RCP<OSMV> B_in_nonconst = getRowMapMultiVector (B, true);
        *B_in_nonconst = B;
        B_in = rcp_const_cast<const OSMV> (B_in_nonconst);

        TPETRA_EFFICIENCY_WARNING
          (! B.isConstantStride (), std::runtime_error,
           "gaussSeidelCopy: The current implementation requires that B have "
           "constant stride.  Since B does not have constant stride, we had to "
           "copy it into a separate constant-stride multivector.  This is a "
           "limitation of the current implementation and not your fault, but we "
           "still report it as an efficiency warning for your information.");
      }

      for (int sweep = 0; sweep < numSweeps; ++sweep) {
        if (! importer.is_null () && sweep > 0) {
          // We already did the first Import for the zeroth sweep above.
          X_colMap->doImport (*X_domainMap, *importer, INSERT);
        }

        // Do local Gauss-Seidel.
        if (direction != Symmetric) {
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     localDirection);
        }
        else { // direction == Symmetric
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     Forward);
          // Communicate again before the Backward sweep, if necessary.
          if (! importer.is_null ()) {
            X_colMap->doImport (*X_domainMap, *importer, INSERT);
          }
          matrix_->template localGaussSeidel<OS,OS> (*B_in, *X_colMap, D,
                                                     dampingFactor,
                                                     Backward);
        }
      }

      if (copyBackOutput) {
        deep_copy (X, *X_domainMap); // Copy result back into X.
      }
    }

    /// \brief Whether this Operator's apply() method can apply the
    ///   transpose or conjugate transpose.
    ///
    /// This is always true, since it is true for the CrsMatrix that
    /// this object wraps.
    bool hasTransposeApply() const {
      return true;
    }

    //! The domain Map of this Operator.
    Teuchos::RCP<const map_type> getDomainMap () const {
      return matrix_->getDomainMap ();
    }

    //! The range Map of this Operator.
    Teuchos::RCP<const map_type> getRangeMap () const {
      return matrix_->getRangeMap ();
    }

    //@}

  protected:
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix object.
    const Teuchos::RCP<const crs_matrix_type> matrix_;

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

    /// \brief Apply the transpose or conjugate transpose of the
    ///   matrix to X_in, producing Y_in.
    void
    applyTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X_in,
                    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y_in,
                    Teuchos::ETransp mode,
                    Scalar alpha,
                    Scalar beta) const
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      using Teuchos::null;

      const int myImageID = Teuchos::rank(*matrix_->getComm());

      const size_t numVectors = X_in.getNumVectors();
      // because of Views, it is difficult to determine if X and Y point to the same data.
      // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
      // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
      Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer = matrix_->getGraph()->getImporter();
      Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = matrix_->getGraph()->getExporter();
      // access X indirectly, in case we need to create temporary storage
      Teuchos::RCP<const MV> X;

      // some parameters for below
      const bool Y_is_replicated = !Y_in.isDistributed(),
        Y_is_overwritten = (beta == ST::zero());
      if (Y_is_replicated && myImageID > 0) {
        beta = ST::zero();
      }

      // currently, cannot multiply from multivector of non-constant stride
      if (X_in.isConstantStride() == false && importer==null) {
        // generate a strided copy of X_in
        X = Teuchos::rcp(new MV(X_in));
      }
      else {
        // just temporary, so this non-owning RCP is okay
        X = Teuchos::rcp(&X_in, false);
      }

      // set up import/export temporary multivectors
      if (importer != null) {
        if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
        if (importMV_ == null) {
          importMV_ = Teuchos::rcp( new MV(matrix_->getColMap(),numVectors) );
        }
      }
      if (exporter != null) {
        if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
        if (exportMV_ == null) {
          exportMV_ = Teuchos::rcp( new MV(matrix_->getRowMap(),numVectors) );
        }
      }

      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter != null) {
        {
          exportMV_->doImport(X_in,*exporter,INSERT);
        }
        // multiply out of exportMV_
        X = exportMV_;
      }


      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (importer != null) {
        // Do actual computation
        matrix_->template localMultiply<Scalar, Scalar>(*X, *importMV_, mode, alpha, ST::zero());
        if (Y_is_overwritten) Y_in.putScalar(ST::zero());
        else                  Y_in.scale(beta);
        //
        {
          Y_in.doExport(*importMV_,*importer,ADD);
        }
      }
      // otherwise, multiply into Y
      else {
        // can't multiply in-situ; can't multiply into non-strided multivector
        if (Y_in.isConstantStride() == false || X.getRawPtr() == &Y_in) {
          // generate a strided copy of Y
          MV Y(Y_in);
          matrix_->template localMultiply<Scalar, Scalar>(*X, Y, mode, alpha, beta);
          deep_copy (Y_in, Y);
        }
        else {
          matrix_->template localMultiply<Scalar, Scalar>(*X, Y_in, mode, alpha, beta);
        }
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y_is_replicated) {
        Y_in.reduce();
      }
    }

    //! Apply the matrix (not its transpose) to X_in, producing Y_in.
    void
    applyNonTranspose (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X_in,
                       MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y_in,
                       Scalar alpha,
                       Scalar beta) const
    {
      using Tpetra::Details::ProfilingRegion;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_const_cast;
      using Teuchos::rcpFromRef;
      typedef Export<LocalOrdinal,GlobalOrdinal,Node> export_type;
      typedef Import<LocalOrdinal,GlobalOrdinal,Node> import_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;

      if (alpha == STS::zero ()) {
        if (beta == STS::zero ()) {
          Y_in.putScalar (STS::zero ());
        } else if (beta != STS::one ()) {
          Y_in.scale (beta);
        }
        return;
      }

      // It's possible that X is a view of Y or vice versa.  We don't
      // allow this (apply() requires that X and Y not alias one
      // another), but it's helpful to detect and work around this
      // case.  We don't try to to detect the more subtle cases (e.g.,
      // one is a subview of the other, but their initial pointers
      // differ).  We only need to do this if this matrix's Import is
      // trivial; otherwise, we don't actually apply the operator from
      // X into Y.

      RCP<const import_type> importer = matrix_->getGraph()->getImporter();
      RCP<const export_type> exporter = matrix_->getGraph()->getExporter();

      // If beta == 0, then the output MV will be overwritten; none of
      // its entries should be read.  (Sparse BLAS semantics say that we
      // must ignore any Inf or NaN entries in Y_in, if beta is zero.)
      // This matters if we need to do an Export operation; see below.
      const bool Y_is_overwritten = (beta == STS::zero());

      // We treat the case of a replicated MV output specially.
      const bool Y_is_replicated =
        (! Y_in.isDistributed () && matrix_->getComm ()->getSize () != 1);

      // This is part of the special case for replicated MV output.
      // We'll let each process do its thing, but do an all-reduce at
      // the end to sum up the results.  Setting beta=0 on all
      // processes but Proc 0 makes the math work out for the
      // all-reduce.  (This assumes that the replicated data is
      // correctly replicated, so that the data are the same on all
      // processes.)
      if (Y_is_replicated && matrix_->getComm ()->getRank () > 0) {
        beta = STS::zero();
      }

      // Temporary MV for Import operation.  After the block of code
      // below, this will be an (Imported if necessary) column Map MV
      // ready to give to localMultiply().
      RCP<const MV> X_colMap;
      if (importer.is_null ()) {
        if (! X_in.isConstantStride ()) {
          // Not all sparse mat-vec kernels can handle an input MV with
          // nonconstant stride correctly, so we have to copy it in that
          // case into a constant stride MV.  To make a constant stride
          // copy of X_in, we force creation of the column (== domain)
          // Map MV (if it hasn't already been created, else fetch the
          // cached copy).  This avoids creating a new MV each time.
          RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in, true);
          Tpetra::deep_copy (*X_colMapNonConst, X_in);
          X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
        }
        else {
          // The domain and column Maps are the same, so do the local
          // multiply using the domain Map input MV X_in.
          X_colMap = rcpFromRef (X_in);
        }
      }
      else { // need to Import source (multi)vector
        ProfilingRegion regionImport ("Tpetra::CrsMatrixMultiplyOp::apply: Import");

        // We're doing an Import anyway, which will copy the relevant
        // elements of the domain Map MV X_in into a separate column Map
        // MV.  Thus, we don't have to worry whether X_in is constant
        // stride.
        RCP<MV> X_colMapNonConst = getColumnMapMultiVector (X_in);

        // Import from the domain Map MV to the column Map MV.
        X_colMapNonConst->doImport (X_in, *importer, INSERT);
        X_colMap = rcp_const_cast<const MV> (X_colMapNonConst);
      }

      // Temporary MV for doExport (if needed), or for copying a
      // nonconstant stride output MV into a constant stride MV.  This
      // is null if we don't need the temporary MV, that is, if the
      // Export is trivial (null).
      RCP<MV> Y_rowMap = getRowMapMultiVector (Y_in);

      // If we have a nontrivial Export object, we must perform an
      // Export.  In that case, the local multiply result will go into
      // the row Map multivector.  We don't have to make a
      // constant-stride version of Y_in in this case, because we had to
      // make a constant stride Y_rowMap MV and do an Export anyway.
      if (! exporter.is_null ()) {
        matrix_->template localMultiply<Scalar, Scalar> (*X_colMap, *Y_rowMap,
                                                         Teuchos::NO_TRANS,
                                                         alpha, STS::zero());
        {
          ProfilingRegion regionExport ("Tpetra::CrsMatrixMultiplyOp::apply: Export");

          // If we're overwriting the output MV Y_in completely (beta
          // == 0), then make sure that it is filled with zeros before
          // we do the Export.  Otherwise, the ADD combine mode will
          // use data in Y_in, which is supposed to be zero.
          if (Y_is_overwritten) {
            Y_in.putScalar (STS::zero());
          }
          else {
            // Scale the output MV by beta, so that the Export sums in
            // the mat-vec contribution: Y_in = beta*Y_in + alpha*A*X_in.
            Y_in.scale (beta);
          }
          // Do the Export operation.
          Y_in.doExport (*Y_rowMap, *exporter, ADD);
        }
      }
      else { // Don't do an Export: row Map and range Map are the same.
        //
        // If Y_in does not have constant stride, or if the column Map
        // MV aliases Y_in, then we can't let the kernel write directly
        // to Y_in.  Instead, we have to use the cached row (== range)
        // Map MV as temporary storage.
        //
        // FIXME (mfh 05 Jun 2014, mfh 07 Dec 2018) This test for
        // aliasing only tests if the user passed in the same
        // MultiVector for both X and Y.  It won't detect whether one
        // MultiVector views the other.  We should also check the
        // MultiVectors' raw data pointers.
        if (! Y_in.isConstantStride () || X_colMap.getRawPtr () == &Y_in) {
          // Force creating the MV if it hasn't been created already.
          // This will reuse a previously created cached MV.
          Y_rowMap = getRowMapMultiVector (Y_in, true);

          // If beta == 0, we don't need to copy Y_in into Y_rowMap,
          // since we're overwriting it anyway.
          if (beta != STS::zero ()) {
            Tpetra::deep_copy (*Y_rowMap, Y_in);
          }
          matrix_->template localMultiply<Scalar, Scalar> (*X_colMap,
                                                           *Y_rowMap,
                                                           Teuchos::NO_TRANS,
                                                           alpha, beta);
          Tpetra::deep_copy (Y_in, *Y_rowMap);
        }
        else {
          matrix_->template localMultiply<Scalar, Scalar> (*X_colMap, Y_in,
                                                           Teuchos::NO_TRANS,
                                                           alpha, beta);
        }
      }

      // If the range Map is a locally replicated Map, sum up
      // contributions from each process.  We set beta = 0 on all
      // processes but Proc 0 initially, so this will handle the scaling
      // factor beta correctly.
      if (Y_is_replicated) {
        ProfilingRegion regionReduce ("Tpetra::CrsMatrixMultiplyOp::apply: Reduce Y");
        Y_in.reduce ();
      }
    }

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
                             const bool force = false) const
    {
      using Teuchos::null;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef Import<LocalOrdinal,GlobalOrdinal,Node> import_type;
      typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

      const size_t numVecs = X_domainMap.getNumVectors ();
      RCP<const import_type> importer = matrix_->getGraph ()->getImporter ();
      RCP<const map_type> colMap = matrix_->getColMap ();

      RCP<MV> X_colMap; // null by default

      // If the Import object is trivial (null), then we don't need a
      // separate column Map multivector.  Just return null in that
      // case.  The caller is responsible for knowing not to use the
      // returned null pointer.
      //
      // If the Import is nontrivial, then we do need a separate
      // column Map multivector for the Import operation.  Check in
      // that case if we have to (re)create the column Map
      // multivector.
      if (! importer.is_null () || force) {
        if (importMV_.is_null () || importMV_->getNumVectors () != numVecs) {
          X_colMap = rcp (new MV (colMap, numVecs));

          // Cache the newly created multivector for later reuse.
          importMV_ = X_colMap;
        }
        else { // Yay, we can reuse the cached multivector!
          X_colMap = importMV_;
          // mfh 09 Jan 2013: We don't have to fill with zeros first,
          // because the Import uses INSERT combine mode, which overwrites
          // existing entries.
          //
          //X_colMap->putScalar (STS::zero ());
        }
      }
      return X_colMap;
    }

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
                          const bool force = false) const
    {
      using Teuchos::null;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef Export<LocalOrdinal,GlobalOrdinal,Node> export_type;
      typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

      const size_t numVecs = Y_rangeMap.getNumVectors ();
      RCP<const export_type> exporter = matrix_->getGraph ()->getExporter ();
      RCP<const map_type> rowMap = matrix_->getRowMap ();

      RCP<MV> Y_rowMap; // null by default

      // If the Export object is trivial (null), then we don't need a
      // separate row Map multivector.  Just return null in that case.
      // The caller is responsible for knowing not to use the returned
      // null pointer.
      //
      // If the Export is nontrivial, then we do need a separate row
      // Map multivector for the Export operation.  Check in that case
      // if we have to (re)create the row Map multivector.
      if (! exporter.is_null () || force) {
        if (exportMV_.is_null () || exportMV_->getNumVectors () != numVecs) {
          Y_rowMap = rcp (new MV (rowMap, numVecs));
          exportMV_ = Y_rowMap; // Cache the newly created MV for later reuse
        }
        else { // Yay, we can reuse the cached multivector!
          Y_rowMap = exportMV_;
        }
      }
      return Y_rowMap;
    }
  };

  /// \brief Non-member function to create a CrsMatrixMultiplyOp.
  /// \relatesalso CrsMatrixMultiplyOp
  ///
  /// The function has the same template parameters of CrsMatrixMultiplyOp.
  ///
  /// \param A [in] The CrsMatrix instance to wrap in an CrsMatrixMultiplyOp.
  /// \return The CrsMatrixMultiplyOp wrapper for the given CrsMatrix.
  template <class OpScalar,
            class MatScalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  Teuchos::RCP<
    CrsMatrixMultiplyOp<OpScalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node> >
  createCrsMatrixMultiplyOp (const Teuchos::RCP<
    const CrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
  {
    typedef CrsMatrixMultiplyOp<OpScalar, MatScalar, LocalOrdinal,
      GlobalOrdinal, Node> op_type;
    return Teuchos::rcp (new op_type (A));
  }

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXMULTIPLYOP_HPP

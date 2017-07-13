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

#ifndef TPETRA_CRSMATRIXSOLVEOP_HPP
#define TPETRA_CRSMATRIXSOLVEOP_HPP

/// \file Tpetra_CrsMatrixSolveOp.hpp
///
/// Declaration and definition of Tpetra::CrsMatrixSolveOp and its
/// nonmember constructor Tpetra::createCrsMatrixSolveOp.

#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra {

  /// \class CrsMatrixSolveOp
  /// \brief Wrap a CrsMatrix instance's triangular solve in an Operator.
  ///
  /// \tparam Scalar Same as the first template parameter of Operator.
  ///   The type of the entries of the MultiVector input and output of
  ///   apply().  Not necessarily the same as the first template
  ///   parameter of the CrsMatrix used to create this object.
  /// \tparam MatScalar Same as the first template parameter of
  ///   CrsMatrix.  The type of the entries of the sparse matrix.  Not
  ///   necessarily the same as the type of the entries of the
  ///   MultiVector input and output of apply().
  /// \tparam LocalOrdinal Same as the second template parameter of
  ///   CrsMatrix and Operator.
  /// \tparam GlobalOrdinal Same as the third template parameter of
  ///   CrsMatrix and Operator.
  /// \tparam Node Same as the fourth template parameter of CrsMatrix
  ///   and Operator.
  ///
  /// This class' apply() method does a "local" triangular solve.
  /// "Local" is in quotes because apply() does the same communication
  /// (Import and Export) operations that CrsMatrix's apply() method
  /// would do for a sparse matrix-vector multiply, but the triangular
  /// solve is restricted to each process' part of the data.  Thus, it
  /// is not a triangular solve of a fully distributed triangular
  /// matrix.
  ///
  /// Here are some situations where this operation is useful:
  /// - Your sparse matrix A only lives in one MPI process, and you
  ///   have a factorization of it (either complete or incomplete).
  /// - Domain decomposition, where each MPI process owns one subdomain
  /// - Coarse-grid solves in algebraic multigrid
  /// - Mixed-precision operations, where the type <tt>MatScalar</tt>
  ///   of entries in the matrix differs from the type <tt>Scalar</tt>
  ///   of entries in the MultiVector input and output of apply().
  template <class Scalar,
            class MatScalar = Scalar,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
  class CrsMatrixSolveOp :
    public Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    //! The specialization of CrsMatrix which this class wraps.
    typedef CrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    //! The specialization of Map which this class uses.
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    //! @name Constructor and destructor
    //@{

    //! Constructor; takes a CrsMatrix to use for local triangular solves.
    CrsMatrixSolveOp (const Teuchos::RCP<const crs_matrix_type>& A) :
      matrix_ (A)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~CrsMatrixSolveOp () {}

    //@}
    //! @name Implementation of Operator
    //@{

    /// \brief Compute \f$Y = \beta Y + \alpha B X\f$, where \f$B X\f$
    ///   represents the result of the local triangular solve.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const
    {
      typedef Teuchos::ScalarTraits<Scalar> STOS;
      const char prefix[] = "Tpetra::CrsMatrixSolveOp::apply: ";

      TEUCHOS_TEST_FOR_EXCEPTION
        (! matrix_->isFillComplete (), std::runtime_error,
         prefix << "Underlying matrix is not fill complete.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (! matrix_->isLowerTriangular () && ! matrix_->isUpperTriangular (),
         std::runtime_error, prefix << "The matrix is neither lower nor upper "
         "triangular.  Remember that in Tpetra, triangular-ness is a local "
         "(per MPI process) property.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
         prefix << "X and Y must have the same number of columns (vectors).  "
         "X.getNumVectors() = " << X.getNumVectors ()
         << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (alpha != STOS::one () || beta != STOS::zero (), std::logic_error,
         prefix << "The case alpha != 1 or beta != 0 has not yet been "
         "implemented.  Please speak with the Tpetra developers.");
      if (mode == Teuchos::NO_TRANS) {
        applyNonTranspose (X,Y);
      } else if (mode == Teuchos::TRANS || mode == Teuchos::CONJ_TRANS) {
        applyTranspose (X, Y, mode);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, prefix << "The 'mode' argument has an "
           "invalid value " << mode << ".  Valid values are Teuchos::NO_TRANS="
           << Teuchos::NO_TRANS << ", Teuchos::TRANS=" << Teuchos::TRANS << ", "
           "and Teuchos::CONJ_TRANS=" << Teuchos::CONJ_TRANS << ".");
      }
    }

    //! Whether apply() can solve with the (conjugate) transpose of the matrix.
    bool hasTransposeApply () const {
      return true;
    }

    /// \brief The domain Map of this operator.
    /// This is the range map of the underlying CrsMatrix.
    Teuchos::RCP<const map_type> getDomainMap () const
    {
      return matrix_->getRangeMap ();
    }

    /// \brief The range Map of this operator.
    /// This is the domain Map of the underlying CrsMatrix.
    Teuchos::RCP<const map_type> getRangeMap () const {
      return matrix_->getDomainMap ();
    }

    //@}
  protected:
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix.
    const Teuchos::RCP<const crs_matrix_type> matrix_;

    //! Cached temporary destination of Import operation in apply().
    mutable Teuchos::RCP<MV> importMV_;
    //! Cached temporary source of Export operation in apply().
    mutable Teuchos::RCP<MV> exportMV_;

    //! Do the non-transpose solve.
    void applyNonTranspose (const MV& X_in, MV& Y_in) const
    {
      using Teuchos::NO_TRANS;
      using Teuchos::null;
      typedef Teuchos::ScalarTraits<Scalar> ST;

      // Solve U X = Y  or  L X = Y
      // X belongs to domain map, while Y belongs to range map

      const size_t numVectors = X_in.getNumVectors();
      Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
        matrix_->getGraph ()->getImporter ();
      Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter =
        matrix_->getGraph ()->getExporter ();
      Teuchos::RCP<const MV> X;

      // it is okay if X and Y reference the same data, because we can
      // perform a triangular solve in-situ.  however, we require that
      // column access to each is strided.

      // set up import/export temporary multivectors
      if (importer != null) {
        if (importMV_ != null && importMV_->getNumVectors () != numVectors) {
          importMV_ = null;
        }
        if (importMV_ == null) {
          importMV_ = Teuchos::rcp (new MV (matrix_->getColMap (), numVectors));
        }
      }
      if (exporter != null) {
        if (exportMV_ != null && exportMV_->getNumVectors () != numVectors) {
          exportMV_ = null;
        }
        if (exportMV_ == null) {
          exportMV_ = Teuchos::rcp (new MV (matrix_->getRowMap (), numVectors));
        }
      }

      // solve(NO_TRANS): RangeMap -> DomainMap
      // lclMatSolve_: RowMap -> ColMap
      // importer: DomainMap -> ColMap
      // exporter: RowMap -> RangeMap
      //
      // solve = reverse(exporter)  o   lclMatSolve_  o reverse(importer)
      //         RangeMap   ->    RowMap     ->     ColMap         ->    DomainMap
      //
      // If we have a non-trivial exporter, we must import elements that
      // are permuted or are on other processors
      if (exporter != null) {
        exportMV_->doImport (X_in, *exporter, INSERT);
        X = exportMV_;
      }
      else if (! X_in.isConstantStride ()) {
        // cannot handle non-constant stride right now
        // generate a copy of X_in
        X = Teuchos::rcp (new MV (X_in));
      }
      else {
        // just temporary, so this non-owning RCP is okay
        X = Teuchos::rcpFromRef (X_in);
      }

      // If we have a non-trivial importer, we must export elements that
      // are permuted or belong to other processes.  We will compute
      // solution into the to-be-exported MV.
      if (importer != null) {
        matrix_->template localSolve<Scalar, Scalar> (*X, *importMV_, NO_TRANS);
        // Make sure target is zero: necessary because we are adding.
        Y_in.putScalar (ST::zero ());
        Y_in.doExport (*importMV_, *importer, ADD);
      }
      // otherwise, solve into Y
      else {
        // can't solve into non-strided multivector
        if (! Y_in.isConstantStride ()) {
          // generate a strided copy of Y
          MV Y (Y_in);
          matrix_->template localSolve<Scalar, Scalar> (*X, Y, NO_TRANS);
          Tpetra::deep_copy (Y_in, Y);
        }
        else {
          matrix_->template localSolve<Scalar, Scalar> (*X, Y_in, NO_TRANS);
        }
      }
    }

    //! Do the transpose or conjugate transpose solve.
    void applyTranspose (const MV& X_in, MV& Y_in, const Teuchos::ETransp mode) const
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      using Teuchos::null;

      TEUCHOS_TEST_FOR_EXCEPTION
        (mode != Teuchos::TRANS && mode != Teuchos::CONJ_TRANS, std::logic_error,
         "Tpetra::CrsMatrixSolveOp::applyTranspose: mode is neither TRANS nor "
         "CONJ_TRANS.  Should never get here!  Please report this bug to the "
         "Tpetra developers.");

      const size_t numVectors = X_in.getNumVectors();
      Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
        matrix_->getGraph ()->getImporter ();
      Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter =
        matrix_->getGraph ()->getExporter ();
      Teuchos::RCP<const MV> X;

      // it is okay if X and Y reference the same data, because we can
      // perform a triangular solve in-situ.  however, we require that
      // column access to each is strided.

      // set up import/export temporary multivectors
      if (importer != null) {
        if (importMV_ != null && importMV_->getNumVectors() != numVectors) {
          importMV_ = null;
        }
        if (importMV_ == null) {
          importMV_ = Teuchos::rcp( new MV(matrix_->getColMap(),numVectors) );
        }
      }
      if (exporter != null) {
        if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) {
          exportMV_ = null;
        }
        if (exportMV_ == null) {
          exportMV_ = Teuchos::rcp( new MV(matrix_->getRowMap(),numVectors) );
        }
      }

      // solve(TRANS): DomainMap -> RangeMap
      // lclMatSolve_(TRANS): ColMap -> RowMap
      // importer: DomainMap -> ColMap
      // exporter: RowMap -> RangeMap
      //
      // solve = importer o   lclMatSolve_  o  exporter
      //         Domainmap -> ColMap     ->      RowMap -> RangeMap
      //
      // If we have a non-trivial importer, we must import elements that
      // are permuted or are on other processes.
      if (importer != null) {
        importMV_->doImport(X_in,*importer,INSERT);
        X = importMV_;
      }
      else if (X_in.isConstantStride() == false) {
        // cannot handle non-constant stride right now
        // generate a copy of X_in
        X = Teuchos::rcp(new MV(X_in));
      }
      else {
        // just temporary, so this non-owning RCP is okay
        X = Teuchos::rcpFromRef (X_in);
      }


      // If we have a non-trivial exporter, we must export elements that
      // are permuted or belong to other processes.  We will compute
      // solution into the to-be-exported MV; get a view.
      if (exporter != null) {
        matrix_->template localSolve<Scalar, Scalar> (*X, *exportMV_,
                                                      Teuchos::CONJ_TRANS);
        // Make sure target is zero: necessary because we are adding
        Y_in.putScalar(ST::zero());
        Y_in.doExport(*importMV_, *importer, ADD);
      }
      // otherwise, solve into Y
      else {
        if (Y_in.isConstantStride() == false) {
          // generate a strided copy of Y
          MV Y(Y_in);
          matrix_->template localSolve<Scalar, Scalar> (*X, Y, Teuchos::CONJ_TRANS);
          Y_in = Y;
        }
        else {
          matrix_->template localSolve<Scalar, Scalar> (*X, Y_in, Teuchos::CONJ_TRANS);
        }
      }
    }
  };

  /// \brief Nonmember function that wraps a CrsMatrix in a CrsMatrixSolveOp.
  /// \relatesalso CrsMatrixSolveOp
  ///
  /// The function has the same template parameters of CrsMatrixSolveOp.
  ///
  /// \param A [in] The CrsMatrix instance to wrap in an CrsMatrixSolveOp.
  /// \return The CrsMatrixSolveOp wrapper for the given CrsMatrix.
  template<class OpScalar,
           class MatScalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  Teuchos::RCP<CrsMatrixSolveOp<OpScalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node> >
  createCrsMatrixSolveOp (const Teuchos::RCP<const CrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
  {
    return Teuchos::rcp (new CrsMatrixSolveOp<OpScalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node> (A));
  }

} // namespace Tpetra

#endif // TPETRA_CRSMATRIXSOLVEOP_HPP

// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_HPP

/// \file Tpetra_CrsMatrixMultiplyOp.hpp
///
/// Declaration and definition of Tpetra::CrsMatrixMultiplyOp and its
/// nonmember constructor Tpetra::createCrsMatrixMultiplyOp.

#include "Tpetra_CrsMatrixMultiplyOp_fwd.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_LocalCrsMatrixOperator.hpp"

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
  /// The resulting output (Multi)Vector will be computed at its own
  /// precision.
  ///
  /// \tparam Scalar The type of the entries of the input and output
  ///   MultiVector (see apply()).  Same as the first template
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
  ///
  /// If you encounter link errors when Scalar != MatScalar, try
  /// including <tt>Tpetra_LocalCrsMatrixOperator_def.hpp</tt> after
  /// including this header file.
  template <class Scalar,
            class MatScalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class CrsMatrixMultiplyOp :
    public Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    //! The specialization of CrsMatrix which this class wraps.
    using crs_matrix_type =
      CrsMatrix<MatScalar, LocalOrdinal, GlobalOrdinal, Node>;
    //! The specialization of Map which this class uses.
    using map_type = Map<LocalOrdinal, GlobalOrdinal, Node>;

  private:
    using local_matrix_device_type =
      typename crs_matrix_type::local_matrix_device_type;

  public:
    //! @name Constructor and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param A [in] The CrsMatrix to wrap as an
    ///   <tt>Operator<Scalar, ...></tt>.
    CrsMatrixMultiplyOp (const Teuchos::RCP<const crs_matrix_type>& A) :
      matrix_ (A),
      localMultiply_ (std::make_shared<local_matrix_device_type> (
                                       A->getLocalMatrixDevice ()))
    {}

    //! Destructor (virtual for memory safety of derived classes).
    ~CrsMatrixMultiplyOp () override = default;

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute <tt>Y = beta*Y + alpha*Op(A)*X</tt>, where
    ///   <tt>Op(A)</tt> is either A, \f$A^T\f$, or \f$A^H\f$.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const override
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

    /// \brief Whether this Operator's apply() method can apply the
    ///   transpose or conjugate transpose.
    ///
    /// This is always true, since it is true for the CrsMatrix that
    /// this object wraps.
    bool hasTransposeApply() const override {
      return true;
    }

    //! The domain Map of this Operator.
    Teuchos::RCP<const map_type> getDomainMap () const override {
      return matrix_->getDomainMap ();
    }

    //! The range Map of this Operator.
    Teuchos::RCP<const map_type> getRangeMap () const override {
      return matrix_->getRangeMap ();
    }

    //@}

  protected:
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix object.
    const Teuchos::RCP<const crs_matrix_type> matrix_;

    //! Implementation of local sparse matrix-vector multiply.
    LocalCrsMatrixOperator<Scalar, MatScalar, typename
                           crs_matrix_type::device_type> localMultiply_;

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
      using Teuchos::null;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
      using import_type = Import<LocalOrdinal, GlobalOrdinal, Node>;
      using STS = Teuchos::ScalarTraits<Scalar>;
      typedef typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>:: dual_view_type::t_dev nonconst_view_type;

      const size_t numVectors = X_in.getNumVectors();
      // because of Views, it is difficult to determine if X and Y point to the same data.
      // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
      // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
      RCP<const import_type> importer = matrix_->getGraph()->getImporter();
      RCP<const export_type> exporter = matrix_->getGraph()->getExporter();

      // some parameters for below
      const bool Y_is_replicated = ! Y_in.isDistributed ();
      const bool Y_is_overwritten = (beta == STS::zero ());
      const int myRank = matrix_->getComm ()->getRank ();
      if (Y_is_replicated && myRank != 0) {
        beta = STS::zero ();
      }

      // access X indirectly, in case we need to create temporary storage
      RCP<const MV> X;
      // currently, cannot multiply from multivector of non-constant stride
      if (! X_in.isConstantStride () && importer == null) {
        // generate a strided copy of X_in
        X = Teuchos::rcp (new MV (X_in, Teuchos::Copy));
      }
      else {
        // just temporary, so this non-owning RCP is okay
        X = Teuchos::rcpFromRef (X_in);
      }

      // set up import/export temporary multivectors
      if (importer != null) {
        if (importMV_ != null && importMV_->getNumVectors () != numVectors) {
          importMV_ = null;
        }
        if (importMV_ == null) {
          importMV_ = rcp (new MV (matrix_->getColMap (), numVectors));
        }
      }
      if (exporter != null) {
        if (exportMV_ != null && exportMV_->getNumVectors () != numVectors) {
          exportMV_ = null;
        }
        if (exportMV_ == null) {
          exportMV_ = rcp (new MV (matrix_->getRowMap (), numVectors));
        }
      }

      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter != null) {
        exportMV_->doImport(X_in,*exporter,INSERT);
        X = exportMV_; // multiply out of exportMV_
      }

      auto X_lcl = X->getLocalViewDevice(Access::ReadOnly);

      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (importer != null) {
        // Beta is zero here, so we clobber Y_lcl
        auto Y_lcl = importMV_->getLocalViewDevice(Access::OverwriteAll);

        localMultiply_.apply (X_lcl, Y_lcl, mode, alpha, STS::zero ());
        if (Y_is_overwritten) {
          Y_in.putScalar (STS::zero ());
        }
        else {
          Y_in.scale (beta);
        }
        Y_in.doExport (*importMV_, *importer, ADD_ASSIGN);
      }
      // otherwise, multiply into Y
      else {
        // can't multiply in-situ; can't multiply into non-strided multivector
        if (! Y_in.isConstantStride () || X.getRawPtr () == &Y_in) {
          // generate a strided copy of Y
          MV Y (Y_in, Teuchos::Copy);
          nonconst_view_type Y_lcl;
          if(Y_is_overwritten) Y_lcl = Y.getLocalViewDevice(Access::OverwriteAll);
          else                 Y_lcl = Y.getLocalViewDevice(Access::ReadWrite); 

          localMultiply_.apply (X_lcl, Y_lcl, mode, alpha, beta);
          Tpetra::deep_copy (Y_in, Y);
        }
        else {
          nonconst_view_type Y_lcl;
          if(Y_is_overwritten) Y_lcl = Y_in.getLocalViewDevice(Access::OverwriteAll);
          else                 Y_lcl = Y_in.getLocalViewDevice(Access::ReadWrite); 

          localMultiply_.apply (X_lcl, Y_lcl, mode, alpha, beta);
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
      using Teuchos::NO_TRANS;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_const_cast;
      using Teuchos::rcpFromRef;
      typedef Export<LocalOrdinal,GlobalOrdinal,Node> export_type;
      typedef Import<LocalOrdinal,GlobalOrdinal,Node> import_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>:: dual_view_type::t_dev nonconst_view_type;

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
      // ready to give to localMultiply_.apply(...).
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
      
      auto X_lcl = X_colMap->getLocalViewDevice(Access::ReadOnly);

      // If we have a nontrivial Export object, we must perform an
      // Export.  In that case, the local multiply result will go into
      // the row Map multivector.  We don't have to make a
      // constant-stride version of Y_in in this case, because we had to
      // make a constant stride Y_rowMap MV and do an Export anyway.
      if (! exporter.is_null ()) {
        auto Y_lcl = Y_rowMap->getLocalViewDevice(Access::OverwriteAll);

        localMultiply_.apply (X_lcl, Y_lcl, NO_TRANS,
                              alpha, STS::zero ());
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
          Y_in.doExport (*Y_rowMap, *exporter, ADD_ASSIGN);
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
          nonconst_view_type Y_lcl;
          if(Y_is_overwritten) Y_lcl = Y_rowMap->getLocalViewDevice(Access::OverwriteAll);
          else                 Y_lcl = Y_rowMap->getLocalViewDevice(Access::ReadWrite); 

          localMultiply_.apply (X_lcl, Y_lcl, NO_TRANS, alpha, beta);
          Tpetra::deep_copy (Y_in, *Y_rowMap);
        }
        else {
          nonconst_view_type Y_lcl;
          if(Y_is_overwritten) Y_lcl = Y_in.getLocalViewDevice(Access::OverwriteAll);
          else                 Y_lcl = Y_in.getLocalViewDevice(Access::ReadWrite); 

          localMultiply_.apply (X_lcl, Y_lcl, NO_TRANS, alpha, beta);
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

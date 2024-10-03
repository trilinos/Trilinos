// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TPETRA_LOW_RANK_UPDATE_ROW_MATRIX_HPP
#define LOCA_TPETRA_LOW_RANK_UPDATE_ROW_MATRIX_HPP

#include "Tpetra_RowMatrix.hpp" // base class
#include "Tpetra_Operator.hpp" // base class
#include "NOX_TpetraTypedefs.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace LOCA {
  class GlobalData;
}

namespace LOCA {
  namespace Tpetra {

    /*!
     * \brief A Tpetra row matrix for implementing the operator
     * \f$P = J + U V^T\f$.
     */
    /*!
     * This class implements the Tpetra::RowMatrix interface for
     * \f$P = J + U V^T\f$ where \f$J\f$ is an Tpetra::RowMatrix and
     * \f$U\f$ and \f$V\f$ are Tpetra::MultiVectors.  It is derived
     * from LOCA::Tpetra::LowRankUpdateOp to implement the Tpetra::Operator
     * interface.  The interface here implements the Tpetra::RowMatrix
     * interface when the matrix \f$J\f$ is itself a row matrix.  This
     * allows preconditioners to be computed and scaling in linear systems
     * to be performed when using this operator.  The implementation
     * here merely adds the corresponding entries for \f$U V^T\f$ to the
     * rows of \f$J\f$.  Note however this is only an approximation to the
     * true matrix \f$J + U V^T\f$ (which is dense).
     *
     * This class assumes \f$U\f$ and \f$V\f$ have the same distribution
     * as the rows of \f$J\f$.
     */
    class LowRankUpdateRowMatrix :
      virtual public ::Tpetra::RowMatrix<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>,
      virtual public ::Tpetra::Operator<NOX::Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType> {
    public:

      using scalar_type = NOX::Scalar;
      using local_ordinal_type = NOX::LocalOrdinal;
      using global_ordinal_type = NOX::GlobalOrdinal;
      using node_type = NOX::NodeType;
      using mag_type = typename Kokkos::ArithTraits<NOX::Scalar>::mag_type;

      //! Constructor
      /*!
       * \param global_data [in] The global data object
       * \param jacRowMatrix [in] Jacobian operator J as a row matrix
       * \param U_multiVec [in] Multivector representing U
       * \param V_multiVec [in] Multivector representing V
       * \param setup_for_solve [in] Setup data structures for ApplyInverse()
       * \param include_UV_terms [in] Include \f$U V^T\f$ terms in RowMatrix
       *        routines ExtractRowCopy(), ExtactDiagonalCopy(), InvRowSums(),
       *        InvColSums(), NormInf() and NormOne().
       */
      LowRankUpdateRowMatrix(const Teuchos::RCP<LOCA::GlobalData>& global_data,
                             const Teuchos::RCP<NOX::TRowMatrix>& jacRowMatrix,
                             const Teuchos::RCP<NOX::TMultiVector>& U_multiVec,
                             const Teuchos::RCP<NOX::TMultiVector>& V_multiVec,
                             bool setup_for_solve,
                             bool include_UV_terms);

      virtual ~LowRankUpdateRowMatrix() = default;

      //***************************************
      // Derived from Tpetra::RowMatrix interface
      //***************************************
      virtual local_ordinal_type getBlockSize () const override;
      virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const override;
      virtual Teuchos::RCP<const NOX::TMap> getRowMap() const override;
      virtual Teuchos::RCP<const NOX::TMap> getColMap() const override;
      virtual Teuchos::RCP<const NOX::TRowGraph> getGraph() const override;
      virtual ::Tpetra::global_size_t getGlobalNumRows() const override;
      virtual ::Tpetra::global_size_t getGlobalNumCols() const override;
      virtual size_t getLocalNumRows() const override;
      virtual size_t getLocalNumCols() const override;
      virtual NOX::GlobalOrdinal getIndexBase() const override;
      virtual ::Tpetra::global_size_t getGlobalNumEntries() const override;
      virtual size_t getLocalNumEntries() const override;
      virtual size_t getNumEntriesInGlobalRow (NOX::GlobalOrdinal globalRow) const override;
      virtual size_t getNumEntriesInLocalRow (NOX::LocalOrdinal localRow) const override;
      virtual size_t getGlobalMaxNumRowEntries () const override;
      virtual size_t getLocalMaxNumRowEntries () const override;
      virtual bool hasColMap () const override;
      virtual bool isLocallyIndexed() const override;
      virtual bool isGloballyIndexed() const override;
      virtual bool isFillComplete() const override;
      virtual bool supportsRowViews() const override;
      virtual void
      getGlobalRowCopy (NOX::GlobalOrdinal GlobalRow,
                        NOX::TRowMatrix::nonconst_global_inds_host_view_type &Indices,
                        NOX::TRowMatrix::nonconst_values_host_view_type &Values,
                        size_t &NumEntries) const override;
      virtual void
      getLocalRowCopy (NOX::LocalOrdinal LocalRow,
                        NOX::TRowMatrix::nonconst_local_inds_host_view_type &Indices,
                        NOX::TRowMatrix::nonconst_values_host_view_type &Values,
                       size_t &NumEntries) const override;
      virtual void
      getGlobalRowView (NOX::GlobalOrdinal GlobalRow,
                        NOX::TRowMatrix::global_inds_host_view_type &Indices,
                        NOX::TRowMatrix::values_host_view_type &Values) const override;
      virtual void
      getLocalRowView (NOX::LocalOrdinal LocalRow,
                       NOX::TRowMatrix::local_inds_host_view_type &Indices,
                       NOX::TRowMatrix::values_host_view_type &Values) const override;

      // Use the default implementation!
      // virtual NOX::LocalOrdinal
      // getLocalRowViewRaw (const NOX::LocalOrdinal lclRow,
      //                     NOX::LocalOrdinal& numEnt,
      //                     const NOX::LocalOrdinal*& lclColInds,
      //                     const NOX::Scalar*& vals) const;

      virtual void getLocalDiagCopy (NOX::TVector& diag) const override;
      virtual void leftScale (const NOX::TVector& x) override;
      virtual void rightScale (const NOX::TVector& x) override;
      virtual mag_type getFrobeniusNorm() const override;

      // Use the default implementation!
      // virtual Teuchos::RCP<NOX::TRowMatrix>
      // add (const NOX::Scalar& alpha,
      //      const NOX::TRowMatrix& A,
      //      const NOX::Scalar& beta,
      //      const Teuchos::RCP<const NOX::TMap>& domainMap = Teuchos::null,
      //      const Teuchos::RCP<const NOX::TMap>& rangeMap = Teuchos::null,
      //      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

      //***************************************
      // Derived from Tpetra::Operator interface
      //***************************************
      virtual Teuchos::RCP<const NOX::TMap> getDomainMap() const override;
      virtual Teuchos::RCP<const NOX::TMap> getRangeMap() const override;
      virtual void apply(const NOX::TMultiVector &X,
                         NOX::TMultiVector &Y,
                         Teuchos::ETransp mode = Teuchos::NO_TRANS,
                         NOX::Scalar alpha = Teuchos::ScalarTraits<NOX::Scalar>::one(),
                         NOX::Scalar beta = Teuchos::ScalarTraits<NOX::Scalar>::zero()) const override;

    protected:

      //! Compute \c MyRow, \c MyCol entry of \f$U V^T\f$. Views are local Kokkos view types.
      template<typename ViewType>
      KOKKOS_INLINE_FUNCTION
      NOX::Scalar computeUV(int MyRow, int MyCol) const;

    protected:

      //! Stores row matrix representing J
      Teuchos::RCP<NOX::TRowMatrix> J_rowMatrix;

      //! Stores pointer to non-const U
      Teuchos::RCP<NOX::TMultiVector> nonconst_U;

      //! Stores pointer to non-const V
      Teuchos::RCP<NOX::TMultiVector> nonconst_V;

      //! Flag indicating whether to include U*V^T terms
      bool includeUV;

      //! Number of columns in U and V
      int m;

      //! Map for U
      const NOX::TMap& U_map;

      //! Map for V
      const NOX::TMap& V_map;

      //! Row map for J
      const NOX::TMap& row_map;

      //! Temporary workspace
      mutable Teuchos::RCP<NOX::TMultiVector> tmpMat;

      //! Locally replicated map
      Teuchos::RCP<NOX::TMap> local_map;

    }; // class LowRankUpdateRowMatrix

  } // namespace Tpetra

} // namespace Tpetra

#endif // LOCA_TPETRA_LOW_RANK_UPDATE_ROW_MATRIX_HPP

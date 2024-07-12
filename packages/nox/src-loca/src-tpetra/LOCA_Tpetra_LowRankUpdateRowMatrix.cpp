// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TPETRA_LOW_RANK_UPDATE_ROW_MATRIX_DEF_HPP
#define LOCA_TPETRA_LOW_RANK_UPDATE_ROW_MATRIX_DEF_HPP

#include "LOCA_Tpetra_LowRankUpdateRowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_Assert.hpp"

namespace LOCA {
  namespace Tpetra {

    LowRankUpdateRowMatrix::
    LowRankUpdateRowMatrix(const Teuchos::RCP<LOCA::GlobalData>& global_data,
                           const Teuchos::RCP<NOX::TRowMatrix>& jacRowMatrix,
                           const Teuchos::RCP<NOX::TMultiVector>& U_multiVec,
                           const Teuchos::RCP<NOX::TMultiVector>& V_multiVec,
                           bool setup_for_solve,
                           bool include_UV_terms) :
      J_rowMatrix(jacRowMatrix),
      nonconst_U(U_multiVec),
      nonconst_V(V_multiVec),
      includeUV(include_UV_terms),
      m(U_multiVec->getNumVectors()),
      U_map(*U_multiVec->getMap()),
      V_map(*V_multiVec->getMap()),
      row_map(*jacRowMatrix->getRowMap())
    {
      local_map = Teuchos::rcp(new NOX::TMap(U_multiVec->getNumVectors(),0,U_multiVec->getMap()->getComm(),::Tpetra::LocallyReplicated));
    }

    LowRankUpdateRowMatrix::local_ordinal_type 
    LowRankUpdateRowMatrix::getBlockSize () const
    {return J_rowMatrix->getBlockSize();}

    Teuchos::RCP<const Teuchos::Comm<int> >
    LowRankUpdateRowMatrix::getComm() const
    {return J_rowMatrix->getComm();}

    Teuchos::RCP<const NOX::TMap>
    LowRankUpdateRowMatrix::getRowMap() const
    {return J_rowMatrix->getRowMap();}

    Teuchos::RCP<const NOX::TMap>
    LowRankUpdateRowMatrix::getColMap() const
    {return J_rowMatrix->getColMap();}

    Teuchos::RCP<const NOX::TRowGraph>
    LowRankUpdateRowMatrix::getGraph() const
    {return J_rowMatrix->getGraph();}

    ::Tpetra::global_size_t LowRankUpdateRowMatrix::getGlobalNumRows() const
    {return J_rowMatrix->getGlobalNumRows();}

    ::Tpetra::global_size_t LowRankUpdateRowMatrix::getGlobalNumCols() const
    {return J_rowMatrix->getGlobalNumCols();}

    size_t LowRankUpdateRowMatrix::getLocalNumRows() const
    {return J_rowMatrix->getLocalNumRows();}

    size_t LowRankUpdateRowMatrix::getLocalNumCols() const
    {return J_rowMatrix->getLocalNumCols();}

    NOX::GlobalOrdinal LowRankUpdateRowMatrix::getIndexBase() const
    {return J_rowMatrix->getIndexBase();}

    ::Tpetra::global_size_t LowRankUpdateRowMatrix::getGlobalNumEntries() const
    {return J_rowMatrix->getGlobalNumEntries();}

    size_t LowRankUpdateRowMatrix::getLocalNumEntries() const
    {return J_rowMatrix->getLocalNumEntries();}

    size_t LowRankUpdateRowMatrix::getNumEntriesInGlobalRow(NOX::GlobalOrdinal globalRow) const
    {return J_rowMatrix->getNumEntriesInGlobalRow(globalRow);}

    size_t LowRankUpdateRowMatrix::getNumEntriesInLocalRow(NOX::LocalOrdinal localRow) const
    {return J_rowMatrix->getNumEntriesInLocalRow(localRow);}

    size_t LowRankUpdateRowMatrix::getGlobalMaxNumRowEntries() const
    {return J_rowMatrix->getGlobalMaxNumRowEntries();}

    size_t LowRankUpdateRowMatrix::getLocalMaxNumRowEntries() const
    {return J_rowMatrix->getLocalMaxNumRowEntries();}

    bool LowRankUpdateRowMatrix::hasColMap() const
    {return J_rowMatrix->hasColMap();}

    bool LowRankUpdateRowMatrix::isLocallyIndexed() const
    {return J_rowMatrix->isLocallyIndexed();}

    bool LowRankUpdateRowMatrix::isGloballyIndexed() const
    {return J_rowMatrix->isGloballyIndexed();}

    bool LowRankUpdateRowMatrix::isFillComplete() const
    {return J_rowMatrix->isFillComplete();}

    bool LowRankUpdateRowMatrix::supportsRowViews() const
    {return J_rowMatrix->supportsRowViews();}

    void
    LowRankUpdateRowMatrix::getGlobalRowCopy(NOX::GlobalOrdinal GlobalRow,
                                             NOX::TRowMatrix::nonconst_global_inds_host_view_type &Indices,
                                             NOX::TRowMatrix::nonconst_values_host_view_type &Values,
                                             size_t &NumEntries) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getGlobalRowCopy() - NOT implemented yet!");
    }

    void
    LowRankUpdateRowMatrix::getLocalRowCopy (NOX::LocalOrdinal LocalRow,
                                             NOX::TRowMatrix::nonconst_local_inds_host_view_type &Indices,
                                             NOX::TRowMatrix::nonconst_values_host_view_type &Values,
                                             size_t &NumEntries) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getLocalRowCopy() - NOT implemented yet!");
    }

    void
    LowRankUpdateRowMatrix::getGlobalRowView (NOX::GlobalOrdinal GlobalRow,
                                              NOX::TRowMatrix::global_inds_host_view_type &indices,
                                              NOX::TRowMatrix::values_host_view_type &values) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getGlobalRowView() - NOT implemented yet!");
    }

    void
    LowRankUpdateRowMatrix::getLocalRowView(NOX::LocalOrdinal LocalRow,
                                            NOX::TRowMatrix::local_inds_host_view_type &indices,
                                            NOX::TRowMatrix::values_host_view_type &values) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getLocalRowView() - NOT implemented yet!");
    }

    // Use the default implementation!
    // NOX::LocalOrdinal
    // LowRankUpdateRowMatrix::getLocalRowViewRaw(const NOX::LocalOrdinal lclRow,
    //                                            NOX::LocalOrdinal& numEnt,
    //                                            const NOX::LocalOrdinal*& lclColInds,
    //                                            const NOX::Scalar*& vals) const
    // {}

    void LowRankUpdateRowMatrix::getLocalDiagCopy(NOX::TVector& diag) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getLocalDiagCopy() - NOT implemented yet!");
    }

    void LowRankUpdateRowMatrix::leftScale(const NOX::TVector& x)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::leftScale() - NOT implemented yet!");
    }

    void LowRankUpdateRowMatrix::rightScale(const NOX::TVector& x)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::rightScale() - NOT implemented yet!");
    }

    LowRankUpdateRowMatrix::mag_type LowRankUpdateRowMatrix::getFrobeniusNorm() const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                                 "ERROR - LOCA::LowRankRowMatrix::getFrobeniusNorm() - NOT implemented yet!");
    }

    // Use the default implementation!
    // Teuchos::RCP<NOX::TRowMatrix>
    // LowRankUpdateRowMatrix::add(const NOX::Scalar& alpha,
    //                              const NOX::TRowMatrix& A,
    //                              const NOX::Scalar& beta,
    //                              const Teuchos::RCP<const NOX::TMap>& domainMap,
    //                              const Teuchos::RCP<const NOX::TMap>& rangeMap,
    //                              const Teuchos::RCP<Teuchos::ParameterList>& params) const
    // {}

    Teuchos::RCP<const NOX::TMap> LowRankUpdateRowMatrix::getDomainMap() const
    {return J_rowMatrix->getDomainMap();}

    Teuchos::RCP<const NOX::TMap> LowRankUpdateRowMatrix::getRangeMap() const
    {return J_rowMatrix->getRangeMap();}

    void LowRankUpdateRowMatrix::apply(const NOX::TMultiVector &X,
                                       NOX::TMultiVector &Y,
                                       Teuchos::ETransp mode,
                                       NOX::Scalar alpha,
                                       NOX::Scalar beta) const
    {
      TEUCHOS_ASSERT(alpha == Teuchos::ScalarTraits<NOX::Scalar>::one());
      TEUCHOS_ASSERT(beta == Teuchos::ScalarTraits<NOX::Scalar>::zero());

      // Number of input vectors
      auto k = X.getNumVectors();

      // Compute J*Input or J^T*input
      J_rowMatrix->apply(X, Y);

      // Create temporary matrix to store V^T*input or U^T*input
      if (tmpMat.is_null() || tmpMat->getNumVectors() != k)
        tmpMat = Teuchos::rcp(new NOX::TMultiVector(local_map, k, false));

      // if (!useTranspose) {
      if (mode == Teuchos::NO_TRANS) {

        // Compute V^T*Input
        // tmpMat->Multiply('T', 'N', 1.0, *V, Input, 0.0);
        tmpMat->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *nonconst_V, X, 0.0);

        // Compute J*Input + U*(V^T*input)
        // Result.Multiply('N', 'N', 1.0, *U, *tmpMat, 1.0);
        Y.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *nonconst_U, *tmpMat, 1.0);

      }
      else {

        // Compute U^T*Input
        // tmpMat->Multiply('T', 'N', 1.0, *U, Input, 0.0);
        tmpMat->multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *nonconst_U, X, 0.0);

        // Compute J^T*Input + V*(U^T*input)
        Y.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *nonconst_V, *tmpMat, 1.0);

      }
    }

    template <typename ViewType>
    KOKKOS_INLINE_FUNCTION
    NOX::Scalar LowRankUpdateRowMatrix::computeUV(int u_row_lid, int v_row_lid) const
    {
      NOX::Scalar val = 0.0;
      auto U_HostView=nonconst_U->getLocalViewHost(::Tpetra::Access::ReadOnly);
      auto V_HostView=nonconst_V->getLocalViewHost(::Tpetra::Access::ReadOnly);

      // val = sum_{k=0}^m U(i,k)*V(j,k)
      for (int k=0; k<m; ++k)
        val += U_HostView(u_row_lid,k) * V_HostView(v_row_lid,k);

      return val;
    }

  } // namespace Tpetra

} // namespace LOCA

#endif // TPETRA_ROWMATRIX_DECL_HPP

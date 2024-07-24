// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_TsqrAdaptor.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include <limits>

namespace { // (anonymous)

  Teuchos::RCP<const Tpetra::Map<>>
  makeMap(const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
  {
    using map_type = Tpetra::Map<>;

    TEUCHOS_ASSERT( ! comm.is_null() );

    const map_type::local_ordinal_type lclNumInds (5);
    const int numProcs = comm->getSize();
    const Tpetra::global_size_t gblNumInds (numProcs * lclNumInds);
    const map_type::global_ordinal_type indexBase (0);
    return Teuchos::rcp(new map_type(gblNumInds, indexBase, comm));
  }

  TEUCHOS_UNIT_TEST( TsqrAdaptor, WhetherItCompiles )
  {
    using MV = Tpetra::MultiVector<>;
    using tsqr_adaptor_type = Tpetra::TsqrAdaptor<MV>;

    out << "TsqrAdaptor type name: "
        << Teuchos::TypeNameTraits<tsqr_adaptor_type>::name()
        << std::endl;
    auto map = makeMap(Tpetra::getDefaultComm());

    const size_t numVecs = 3;
    const bool zeroOut = false;
    MV A(map, numVecs, zeroOut);
    MV Q(map, numVecs);
    A.randomize();

    tsqr_adaptor_type thingie;

    // FIXME (mfh 02 Dec 2019) It's not clear that setParameterList
    // works for an empty but nonnull ParameterList.
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    RCP<ParameterList> tsqrParams;
    thingie.setParameterList(tsqrParams);

    RCP<const ParameterList> validTsqrParams =
      thingie.getValidParameters();
    tsqrParams = rcp(new ParameterList(*validTsqrParams));
    thingie.setParameterList(tsqrParams);

    using LO = MV::local_ordinal_type;
    using SC = MV::scalar_type;
    using dense_matrix_type = Teuchos::SerialDenseMatrix<LO, SC>;
    dense_matrix_type R(numVecs, numVecs, zeroOut);

    for (bool forceNonnegativeDiagonal : {false, true}) {
      thingie.factorExplicit (A, Q, R, forceNonnegativeDiagonal);
      // Random matrices need not necessarily have full rank.
      // Just print the rank so we don't emit a build warning.
      using mag_type = MV::mag_type;
      const mag_type tol = std::numeric_limits<mag_type>::epsilon();
      const int rank = thingie.revealRank(Q, R, tol);
      out << "Current rank: " << rank << std::endl;
    }
  }

} // namespace (anonymous)


// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
//  This test exercises Thyra's TSQR adapter.
//

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include "BelosConfigDefs.hpp"

#  include "Thyra_TsqrAdaptor.hpp"
#  include "Thyra_TpetraThyraWrappers.hpp"

int
main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using map_type = Tpetra::Map<>;
  using mv_type = Tpetra::MultiVector<>;
  using Scalar = typename mv_type::scalar_type;

  bool success = true;

  Tpetra::initialize(&argc, &argv);
  int myRank = 0;
  {
    auto comm = Tpetra::getDefaultComm();
    myRank = comm->getRank();

    // Number of global elements
    const int globalNumRows = 100;
    const int blockSize = 3;
    const int indexBase = 0;

    Teuchos::RCP<const map_type> range_tpetra = Teuchos::rcp(new map_type (globalNumRows, indexBase, comm));
    auto X_tpetra = Teuchos::rcp(new mv_type(range_tpetra, blockSize));

    // Create a Thyra vector space.
    RCP<const Thyra::VectorSpaceBase<Scalar> > range_thyra = Thyra::createVectorSpace<Scalar> (range_tpetra);
    // Create a multivector from the Tpetra::MultiVector.
    RCP<Thyra::MultiVectorBase<Scalar> > X_thyra = Thyra::createMultiVector (X_tpetra, range_thyra);

    // typedef Thyra::TsqrAdaptor<Scalar> tsqr_adapter_type;
  }

  if (success) {
    if (myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  else {
    if (myRank == 0) {
      std::cout << "End Result: TEST FAILED" << std::endl;
    }
    return EXIT_FAILURE;
  }

  Tpetra::finalize();
}

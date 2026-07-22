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

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#endif

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosEpetraAdapter.hpp"

#ifdef HAVE_EPETRA_THYRA
#  include "Thyra_TsqrAdaptor.hpp"
#  include "Thyra_EpetraThyraWrappers.hpp"
#  include "Thyra_EpetraLinearOp.hpp"
#endif // HAVE_EPETRA_THYRA

int
main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  bool success = true;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init (&argc, &argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else // NOT HAVE_MPI
  // If we aren't using MPI, then setup a serial communicator.
  Epetra_SerialComm comm;
#endif // HAVE_MPI
  const int myRank = comm.MyPID ();

  // Number of global elements
  const int globalNumRows = 100;
  const int blockSize = 3;
  const int indexBase = 0;

  RCP<const Epetra_Map> range_epetra (new Epetra_Map (globalNumRows, indexBase, comm));
  RCP<Epetra_MultiVector> X_epetra (new Epetra_MultiVector (*range_epetra, blockSize));

  // // Get update list and number of local equations from newly created Map.
  // int NumMyElements = Map.NumMyElements();
  // std::vector<int> MyGlobalElements(NumMyElements);
  // Map->MyGlobalElements(&MyGlobalElements[0]);

#ifdef HAVE_EPETRA_THYRA
  // Create a Thyra vector space.
  RCP<const Thyra::VectorSpaceBase<double> > range_thyra =
    Thyra::create_VectorSpace (range_epetra);
  // Create a multivector from the Epetra_MultiVector.
  RCP<Thyra::MultiVectorBase<double> > X_thyra =
    Thyra::create_MultiVector (rcp_implicit_cast<Epetra_MultiVector> (X_epetra), range_thyra);

  (void) range_thyra;
  (void) X_thyra;

  typedef Thyra::TsqrAdaptor<double> tsqr_adapter_type;

#endif // HAVE_EPETRA_THYRA

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

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
}

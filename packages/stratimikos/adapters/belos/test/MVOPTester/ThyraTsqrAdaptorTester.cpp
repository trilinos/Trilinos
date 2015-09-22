// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
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

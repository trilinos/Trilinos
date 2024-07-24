// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_checkGlobalError.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TestForException.hpp"
#include <iostream>
#include <stdexcept>

namespace Tpetra {
namespace Details {

void
checkGlobalError(std::ostream& globalOutputStream,
                 const bool localSuccess,
                 const char localErrorMessage[],
                 const char globalErrorMessageHeader[],
                 const Teuchos::Comm<int>& comm)
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;

  int lclGood = localSuccess ? 1 : 0;
  int gblGood = 0;
  reduceAll(comm, REDUCE_MIN, lclGood, outArg(gblGood));
  if (gblGood != 1) {
    const int myRank = comm.getRank();
    if (myRank == 0) {
      globalOutputStream << endl << globalErrorMessageHeader
                         << endl;
    }

    if (localSuccess || localErrorMessage == nullptr) {
      Details::gathervPrint(globalOutputStream, "", comm);
    }
    else {
      std::ostringstream lclMsg;
      lclMsg << endl;
      constexpr int numStars = 60;
      for (int star = 0; star < numStars; ++star) {
        lclMsg << '*';
      }
      lclMsg << endl << "Proc " << myRank << ": "
             << localErrorMessage << endl;
      Details::gathervPrint(globalOutputStream, lclMsg.str(), comm);
    }

#ifdef HAVE_TPETRA_MPI
    (void) MPI_Abort(MPI_COMM_WORLD, -1);
#else
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::runtime_error, "Tpetra reports a global error.");
#endif // HAVE_TPETRA_MPI
  }
}

} // namespace Details
} // namespace Tpetra

/*
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
// ************************************************************************
// @HEADER
*/

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

#ifdef HAVE_TEUCHOS_MPI
    (void) MPI_Abort(MPI_COMM_WORLD, -1);
#else
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::runtime_error, "Tpetra reports a global error.");
#endif // HAVE_TEUCHOS_MPI
  }
}

} // namespace Details
} // namespace Tpetra

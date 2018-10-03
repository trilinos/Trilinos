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
#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;

    // Tpetra's default communicator will use a 1-process comm.
#ifdef HAVE_TPETRA_MPI
  Tpetra::ScopeGuard tpetraScope(&argc, &argv, MPI_COMM_SELF);
#else
  // Not building with MPI, so default comm won't use MPI.
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
#endif // HAVE_TPETRA_MPI

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // With a "serial" communicator, the rank is always 0,
  // and the number of processes is always 1.
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  if (myRank == 0) {
    std::cout << "Total number of processes: " << numProcs << std::endl;
  }

  if (comm->getRank () == 0) {
    // On Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version() << std::endl << std::endl;
  }

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;
}

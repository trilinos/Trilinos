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
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>

int
main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  typedef Tpetra::Map<> map_type;
  typedef typename map_type::global_ordinal_type global_ordinal_type;

  Teuchos::oblackholestream blackHole;
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {

    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();

    // Create sizes for map
    const size_t numLocalEntries = 5;
    const Tpetra::global_size_t numGlobalEntries = comm->getSize() * numLocalEntries;
    const global_ordinal_type indexBase = 0;

    // Create a Map which has the same number of global entries per
    // process as contigMap, but distributes them differently, in
    // round-robin (1-D cyclic) fashion instead of contiguously.
    RCP<const map_type> cyclicMap;
    {
      Array<global_ordinal_type>::size_type numEltsPerProc = 5;
      Array<global_ordinal_type> elementList(numEltsPerProc);
      const int numProcs = comm->getSize();
      const int myRank = comm->getRank();
      for (Array<global_ordinal_type>::size_type k = 0; k < numEltsPerProc; ++k) {
        elementList[k] = myRank + k*numProcs;
      }
      cyclicMap = rcp(new map_type(numGlobalEntries, elementList, indexBase, comm));
    }

    // If there's more than one MPI process in the communicator,
    // then cyclicMap is definitely NOT contiguous.
    TEUCHOS_TEST_FOR_EXCEPTION(
      comm->getSize() > 1 && cyclicMap->isContiguous(),
      std::logic_error,
      "The cyclic Map claims to be contiguous.");

    // Create a contiguous and uniform map to check for sameness and compatibility
    RCP<const map_type> contigMap = rcp(new map_type(numGlobalEntries, indexBase, comm));

    // contigMap and cyclicMap should always be compatible.  However, if
    // the communicator contains more than 1 process, then contigMap and
    // cyclicMap are NOT the same.
    TEUCHOS_TEST_FOR_EXCEPTION(! contigMap->isCompatible(*cyclicMap),
      std::logic_error,
      "contigMap should be compatible with cyclicMap, but it's not.");

    TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() > 1 && contigMap->isSameAs(*cyclicMap),
      std::logic_error,
      "contigMap should be compatible with cyclicMap, but it's not.");

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}

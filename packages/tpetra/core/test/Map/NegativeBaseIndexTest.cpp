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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::Array;
  using Tpetra::global_size_t;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST( Map, Bug5401_NegativeBaseIndex )
  {
    using std::endl;
    using map_type = Tpetra::Map<>;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using size_type = Teuchos::Array<GO>::size_type;

    out << "Bug 5401 (negative index base) test" << endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm ();
    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();

    out << "Comm has " << numProcs << " process"
        << (numProcs != 1 ? "es" : "") << endl;
    TEST_EQUALITY( numProcs, 2 );
    if (numProcs != 2) {
      out << "This test only works when running with exactly "
        "2 MPI processes." << endl;
      return;
    }

    TEST_ASSERT( std::is_signed<GO>::value );
    if (! std::is_signed<GO>::value) {
      out << "This test only works when the default GlobalOrdinal "
        "type is signed." << endl;
      return;
    }

    // failure reading 1x4 matrix under MPI
    const GO numElements = 78;
    const GO baseIndexIsNegOne = -1;
    const global_size_t GINV = Teuchos::OrdinalTraits<global_size_t>::invalid();
    Teuchos::Array<GO> elements (numElements);

    out << "Create array of global indices.  All processes have the same "
        << "global index.  The first global index on all processes is "
        << baseIndexIsNegOne << "." << endl;

    // first global element is -1
    for (size_type i = 0; i < elements.size (); ++i) {
      elements[i] = static_cast<GO> (i - 1);
    }

    //int localMapCtorSuccess = 0;
    RCP<map_type> map (new map_type (GINV, elements (),
                                     baseIndexIsNegOne, comm));
    out << "Process " << myRank << ":" << endl;
    {
      Teuchos::OSTab tab1 (out);
      out << "My number of global indices: " << map->getNodeNumElements ()
          << endl
          << "Global number of global indices: " << map->getGlobalNumElements ()
          << endl
          << "Index base: " << map->getIndexBase () << endl
          << "My min global index: " << map->getMinGlobalIndex () << endl
          << "Global min global index: " << map->getMinAllGlobalIndex () << endl;
    }

    TEST_EQUALITY( map->getNodeNumElements(),
                   static_cast<size_t> (numElements) );
    TEST_EQUALITY( map->getGlobalNumElements(),
                   static_cast<global_size_t> (numElements*numProcs) );
    TEST_EQUALITY( map->getIndexBase(), static_cast<GO> (-1) );
    TEST_EQUALITY( map->getMinGlobalIndex(),    static_cast<GO> (-1) );
    TEST_EQUALITY( map->getMinAllGlobalIndex(), static_cast<GO> (-1) );

    // All procs fail if any proc fails
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

}



// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Teuchos_RCP.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_MapFactory.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

// This file regroups tests that are specific to Xpetra.

namespace {
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Comm;

  bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    using Xpetra::DefaultPlatform;

    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  // Test of correctness for the return value of the function: Map::getRemoteIndexList()
  // (getRemoteIndexList() uses Xpetra::LookupStatus toXpetra(int))
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, getRemoteIndexList, M, LO, GO )
  {
    typedef typename Teuchos::ArrayView<int>::const_iterator IntConstIt;
    typedef typename ArrayRCP<GO>::size_type GOSize;

    RCP<const Comm<int> > comm = getDefaultComm();
    const int    numImages = comm->getSize();
    const size_t numGlobalEntries = numImages*4;

    const M map(numGlobalEntries, 0, comm);

    size_t numID = numGlobalEntries;
    ArrayRCP<GO>  GIDList(numID);
    ArrayRCP<int> nodeIDList(GIDList.size());
    Xpetra::LookupStatus ls;

    // Test 1 (AllIDsPresent)
    for(GOSize i = 0; i < GIDList.size(); i++) { GIDList[i] = i; }
    ls = map.getRemoteIndexList(GIDList(), nodeIDList());
    TEST_EQUALITY(ls, Xpetra::AllIDsPresent);
    out << "nodeIDList(): "; for(IntConstIt it = nodeIDList.begin(); it != nodeIDList.end(); ++it) { out << *it << " "; }

    // Test 2 (IDNotPresent)
    for(GOSize i = 0; i < GIDList.size(); i++) { GIDList[i] = i; }
    GIDList[GIDList.size()-1] = numGlobalEntries+1; // IDNotPresent
    ls = map.getRemoteIndexList(GIDList(), nodeIDList());
    TEST_EQUALITY(ls, Xpetra::IDNotPresent);
    out << "nodeIDList(): "; for(IntConstIt it = nodeIDList.begin(); it != nodeIDList.end(); ++it) { out << *it << " "; }
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL_( M, LO, GO )                        \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, getRemoteIndexList, M, LO, GO )

#ifdef HAVE_XPETRA_TPETRA

#define UNIT_TEST_GROUP_ORDINAL_TPETRA( LO, GO ) \
  typedef Xpetra::TpetraMap<LO,GO> TpetraMap ## LO ## GO;       \
    UNIT_TEST_GROUP_ORDINAL_(TpetraMap ## LO ## GO, LO, GO)

  UNIT_TEST_GROUP_ORDINAL_TPETRA(int , int)

#endif // HAVE_XPETRA_TPETRA

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::EpetraMap EpetraMap;
  UNIT_TEST_GROUP_ORDINAL_(EpetraMap, int , int)
#endif // HAVE_XPETRA_EPETRA

}

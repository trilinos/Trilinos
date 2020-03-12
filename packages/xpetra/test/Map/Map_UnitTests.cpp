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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Teuchos_as.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#include "Tpetra_Details_Behavior.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  bool mapDebugChecksEnabled()
  {
#ifdef HAVE_XPETRA_TPETRA
    return Tpetra::Details::Behavior::debug("Map");
#else
    return false;
#endif // HAVE_XPETRA_TPETRA
  }

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define TEST_IS_COMPATIBLE(m1,m2,is_compat)               \
{                                                         \
    TEST_EQUALITY_CONST(m1.isCompatible(m1), true);       \
    TEST_EQUALITY_CONST(m2.isCompatible(m2), true);       \
    TEST_EQUALITY_CONST(m1.isCompatible(m2), is_compat);  \
    TEST_EQUALITY_CONST(m2.isCompatible(m1), is_compat);  \
}

#define TEST_IS_SAME_AS(m1,m2,is_sameas)               \
{                                                      \
    TEST_EQUALITY_CONST(m1.isSameAs(m1), true);        \
    TEST_EQUALITY_CONST(m2.isSameAs(m2), true);        \
    TEST_EQUALITY_CONST(m1.isSameAs(m2), is_sameas);   \
    TEST_EQUALITY_CONST(m2.isSameAs(m1), is_sameas);   \
}

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    return Teuchos::rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, validConstructor1, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const Xpetra::global_size_t GSTI = 100;

    TEST_NOTHROW(M map(GSTI,0,comm));

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, validConstructor2, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_NOTHROW(M map(numImages, Teuchos::tuple<GO>(myImageID), 0, comm));
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, validConstructor3, M, LO, GO, N )
  {
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA
    // create Kokkos templates
    typedef typename N::device_type device_type;
    typedef typename device_type::execution_space execution_space;
    // create a comm
    auto comm = getDefaultComm();
    const auto numProcs = comm->getSize();
    const auto myRank    = comm->getRank();

    const auto INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    // Uniform and contiguous map
    {
      const int numDofsPerProc = 2;
      const int offset = myRank*numDofsPerProc;

      Kokkos::View<GO*, typename N::device_type> indexList("Xpetra: Map indexList", numDofsPerProc);
      Kokkos::parallel_for("Xpetra: Map unit test",
                           Kokkos::RangePolicy<execution_space>(0, numDofsPerProc),
                           KOKKOS_LAMBDA(const int i) {
                             indexList(i) = GO (offset + i);
                           });
      M m(INVALID, indexList, 0, comm);

      TEST_EQUALITY(m.getGlobalNumElements(),
                    Teuchos::as<Xpetra::global_size_t>(numDofsPerProc*numProcs));

      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getNodeNumElements(), Teuchos::as<size_t>(numDofsPerProc));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMinGlobalIndex(), Teuchos::as<GO>(numDofsPerProc*myRank));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMaxGlobalIndex(), Teuchos::as<GO>(numDofsPerProc*myRank + 1));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }

    // Uniform non-contiguous map
    {
      const int numDofsPerProc = 2;

      Kokkos::View<GO*, typename N::device_type> indexList("Xpetra: Map indexList", numDofsPerProc);
      Kokkos::parallel_for("Xpetra: Map unit test",
                           Kokkos::RangePolicy<execution_space>(0, numDofsPerProc),
                           KOKKOS_LAMBDA(const int i) {
                             indexList(i) = GO (3*numProcs*i + myRank);
                           });
      M m(INVALID, indexList, 0, comm);

      TEST_EQUALITY(m.getGlobalNumElements(),
                    Teuchos::as<Xpetra::global_size_t>(numDofsPerProc*numProcs));

      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getNodeNumElements(), Teuchos::as<size_t>(numDofsPerProc));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMinGlobalIndex(), Teuchos::as<GO>(myRank));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMaxGlobalIndex(), Teuchos::as<GO>(3*numProcs + myRank));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      std::cout << std::endl << "myRank=" << myRank
                << ", getMinGlobalIndex=" << m.getMinGlobalIndex()
                << ", getMaxGlobalIndex=" << m.getMaxGlobalIndex()
                << std::endl;
    }

    // Non-uniform contiguous map
    {
      const int numDofsPerProc = myRank + 1;
      GO offset = 0;
      for(int rankIdx = 0; rankIdx < myRank; ++rankIdx) {
        offset += rankIdx + 1;
      }

      Kokkos::View<GO*, typename N::device_type> indexList("Xpetra: Map indexList", numDofsPerProc);
      Kokkos::parallel_for("Xpetra: Map unit test",
                           Kokkos::RangePolicy<execution_space>(0, numDofsPerProc),
                           KOKKOS_LAMBDA(const int i) {
                             indexList(i) = GO (offset + i);
                           });
      M m(INVALID, indexList, 0, comm);

      Xpetra::global_size_t globalNumElements = 0;
      for(int rankIdx = 0; rankIdx < numProcs; ++rankIdx) {
        globalNumElements += Teuchos::as<Xpetra::global_size_t>(rankIdx + 1);
      }
      TEST_EQUALITY(m.getGlobalNumElements(), globalNumElements);

      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getNodeNumElements(), Teuchos::as<size_t>(numDofsPerProc));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMinGlobalIndex(), Teuchos::as<GO>(offset));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMaxGlobalIndex(), Teuchos::as<GO>(offset + myRank));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }

    // Non-uniform non-contiguous map
    {
      const int numDofsPerProc = myRank + 1;

      Kokkos::View<GO*, typename N::device_type> indexList("Xpetra: Map indexList", numDofsPerProc);
      Kokkos::parallel_for("Xpetra: Map unit test",
                           Kokkos::RangePolicy<execution_space>(0, numDofsPerProc),
                           KOKKOS_LAMBDA(const int i) {
                             indexList(i) = GO (i*numProcs + myRank);
                           });
      M m(INVALID, indexList, 0, comm);

      Xpetra::global_size_t globalNumElements = 0;
      for(int rankIdx = 0; rankIdx < numProcs; ++rankIdx) {
        globalNumElements += Teuchos::as<Xpetra::global_size_t>(rankIdx + 1);
      }
      TEST_EQUALITY(m.getGlobalNumElements(), globalNumElements);

      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getNodeNumElements(), Teuchos::as<size_t>(numDofsPerProc));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMinGlobalIndex(), Teuchos::as<GO>(myRank));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );

      TEST_EQUALITY(m.getMaxGlobalIndex(), Teuchos::as<GO>(myRank*(numProcs + 1)));

      // All procs fail if any proc fails
      globalSuccess_int = -1;
      reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
#endif
#endif
  }

  // This test exercises Tpetra's debug-mode checks.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, invalidConstructor1, M, LO, GO, N )
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using std::endl;

    if (! mapDebugChecksEnabled()) {
      return;
    }
    int lclSuccess = 1;
    int gblSuccess = 0;

    const Xpetra::global_size_t GSTI =
      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    const int totalNumTests = numImages > 1 ? 4 : 1;
    int numPassedTests = 0;

    // bad constructor calls: (num global elements, index base)
    TEST_THROW(M map(GSTI,0,comm), std::invalid_argument);

    lclSuccess = success ? 1 : 0;
    reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Test " << (numPassedTests+1) << " of "
          << totalNumTests << " passed" << endl;
    }
    ++numPassedTests;

    if (numImages == 1) {
      return;
    }

    TEST_THROW(M map((myImageID == 0 ? GSTI : 0),0,comm), std::invalid_argument);

    lclSuccess = success ? 1 : 0;
    reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Test " << (numPassedTests+1) << " of "
          << totalNumTests << " passed" << endl;
    }
    ++numPassedTests;

    TEST_THROW(M map((myImageID == 0 ?  1 : 0),0,comm), std::invalid_argument);

    lclSuccess = success ? 1 : 0;
    reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Test " << (numPassedTests+1) << " of "
          << totalNumTests << " passed" << endl;
    }
    ++numPassedTests;

    TEST_THROW(M map(0,(myImageID == 0 ? 0 : 1), comm), std::invalid_argument);

    lclSuccess = success ? 1 : 0;
    reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Test " << (numPassedTests+1) << " of "
          << totalNumTests << " passed" << endl;
    }
    ++numPassedTests;
  }

  // This test exercises Tpetra's debug-mode checks.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, invalidConstructor2, M, LO, GO, N )
  {
    if (! mapDebugChecksEnabled()) {
      return;
    }

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    // bad constructor calls: (num global elements, num local elements, index base)
    TEST_THROW(M map(1,0,0, comm),  std::invalid_argument);
    if (numImages > 1) {
      TEST_THROW(M map((myImageID == 0 ? GSTI :  1),0,0,comm), std::invalid_argument);
      TEST_THROW(M map((myImageID == 0 ?  1 :  0),0,0,comm), std::invalid_argument);
      TEST_THROW(M map(0,0,(myImageID == 0 ? 0 : 1),comm), std::invalid_argument);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

#if 0 // failing for epetra. Epetra does not throw
#ifdef HAVE_TPETRA_DEBUG
  // This test will only pass in a debug build of Tpetra (HAVE_TPETRA_DEBUG).
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, invalidConstructor3, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    // bad constructor calls: (num global, entry list, index base)
    TEST_THROW(M map(numImages, tuple<GO>(-myImageID), 1, comm), std::invalid_argument); // GID less than iB
    if (numImages > 1) {
      TEST_THROW(M map( 1, tuple<GO>(myImageID+1), 1, comm), std::invalid_argument);    // nG != sum nL
      TEST_THROW(M map((myImageID == 0 ? GSTI :  0),tuple<GO>(myImageID+1),1, comm), std::invalid_argument);
      TEST_THROW(M map(0, tuple<GO>(myImageID+1), (myImageID == 0 ? 0 : 1), comm), std::invalid_argument);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif // HAVE_TPETRA_DEBUG
#endif // end if 0

#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TPETRA_ENABLE_SS_TESTING) && defined(HAVE_MPI)
  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, RogersUnsignedGOBugVerification, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    const int myImageID = comm->getRank();
    const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    Teuchos::RCP<M> m;
    TEST_NOTHROW( m = rcp(new M(GSTI, Teuchos::tuple<size_t>(myImageID), 0, comm)) );
    if (m != Teuchos::null) {
      TEST_EQUALITY( m->getMinAllGlobalIndex(), (size_t)0 );
      TEST_EQUALITY( m->getMaxAllGlobalIndex(), (size_t)numImages-1 );
    }
  }
#endif


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, compatabilityTests, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    // test isCompatible()
    // m1.isCompatible(m2) should be true if m1 and m2 have the same number of global entries and the same number of local entries on
    // corresponding nodes
    // test the following scenarios:
    // * same number of global and local entries on all nodes
    // * same number of global entries, but different number of local entries on every node
    // * same number of global entries, but different number of local entries on some nodes
    // * different number of global entries, different number of local entries
    //
    // for each, also:
    // test symmetry   : m1.isCompatible(m2) <=> m2.isCompatible(m1)
    // test reflexivity: m1.isCompatible(m1), m2.isCompatible(m2)
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, true );
    }
    {
      M m1(GSTI,myImageID+1,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, false);
    }
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      {
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,(myImageID+1)%numImages,0,comm);
        TEST_IS_COMPATIBLE( m1, m2, false);
      }
      if (numImages > 2) {
        // want different num local on a subset of procs
        // image 0 and numImages-1 get map1:numLocal==[0,numImages-1] and map2:numLocal==[numImages-1,0], the others get numLocal==myImageID
        LO mynl1, mynl2;
        if (myImageID == 0) {
          mynl1 = 0;
          mynl2 = numImages-1;
        }
        else if (myImageID == numImages-1) {
          mynl1 = numImages-1;
          mynl2 = 0;
        }
        else {
          mynl1 = mynl2 = myImageID;
        }
        {
          M m1(GSTI,mynl1,0,comm),
            m2(GSTI,mynl2,0,comm);
          TEST_IS_COMPATIBLE( m1, m2, false);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, localMap, M, LO, GO, N )
  {
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
#ifdef HAVE_XPETRA_TPETRA  // Note: get Kokkos interface for Epetra is only available if Tpetra is also enabled!
    // create a comm
    auto comm = getDefaultComm();
    const auto numProcs = comm->getSize();
    const auto myRank    = comm->getRank();

    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();
    if (lib == Xpetra::UseEpetra) {
      out << "Xpetra: localMap only valid when using Tpetra" << std::endl;
      return;
    }

    const auto INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    {
      // Contiguous map
      const int numDofsPerProc = 2;
      const int offset = myRank*numDofsPerProc;

      M m(numDofsPerProc*numProcs, 0/*indexBase*/, comm);
      auto localMap = m.getLocalMap();

      TEST_EQUALITY(localMap.getNodeNumElements(), numDofsPerProc);
      for (int i = 0; i < numDofsPerProc; i++)
        TEST_EQUALITY(localMap.getGlobalElement(i), offset + i);
    }
    {
      // Permuted map
      const int numDofsPerProc = 4;
      const int indexBase = 1;
      const int offset = myRank*numDofsPerProc + indexBase;

      Teuchos::Array<GO> elementList(numDofsPerProc);
      elementList[0] = offset + 0;
      elementList[1] = offset + 1;
      elementList[2] = offset + 3;
      elementList[3] = offset + 2;

      M m(INVALID, elementList, indexBase, comm);
      auto localMap = m.getLocalMap();

      TEST_EQUALITY(localMap.getNodeNumElements(), numDofsPerProc);
      for (int i = 0; i < numDofsPerProc; i++)
        TEST_EQUALITY(localMap.getGlobalElement(i), elementList[i]);
    }
    {
      // Sparse map
      const int numDofsPerProc = 4;
      const int indexBase = 1;
      const int offset = myRank*numDofsPerProc + indexBase;

      Teuchos::Array<GO> elementList(numDofsPerProc);
      elementList[0] = offset + 10;
      elementList[1] = offset + 6;
      elementList[2] = offset + 4;
      elementList[3] = offset + 2;

      M m(INVALID, elementList, indexBase, comm);
      auto localMap = m.getLocalMap();

      TEST_EQUALITY(localMap.getNodeNumElements(), numDofsPerProc);
      for (int i = 0; i < numDofsPerProc; i++)
        TEST_EQUALITY(localMap.getGlobalElement(i), elementList[i]);
    }
#endif
#endif
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, sameasTests, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const Xpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    {
      M m1(GSTI,0,0,comm),
        m2(GSTI,0,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID+1,0,comm);
      TEST_IS_SAME_AS(m1, m2, false);
    }
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      {
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,myImageID+(myImageID==1?1:0),0,comm);
        TEST_IS_SAME_AS(m1, m2, false);
      }
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Map, ContigUniformMap, M, LO, GO, N )
  {
    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with two entries per node
    // this map will have the following entries:
    Teuchos::Array<GO> myGlobal( Teuchos::tuple<GO>(myImageID*2, myImageID*2+1) );
    Teuchos::Array<LO>  myLocal( Teuchos::tuple<LO>(0,1) );

    const size_t numGlobalEntries = numImages*2;
    const GO indexBase = 0;
    M map(numGlobalEntries,indexBase,comm);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getNodeNumElements(), 2);
    TEST_EQUALITY_CONST(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), 1);
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[indexBase]);
    TEST_EQUALITY_CONST(map.getMaxGlobalIndex(), myGlobal[1]);
    TEST_EQUALITY_CONST(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalEntries-1));
    TEST_EQUALITY( map.getLocalElement(myGlobal[0]), myLocal[0] );
    TEST_EQUALITY( map.getLocalElement(myGlobal[1]), myLocal[1] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[0]), myGlobal[0] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[1]), myGlobal[1] );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries), Teuchos::OrdinalTraits<LO>::invalid() );
    TEST_EQUALITY( map.getGlobalElement(2),               Teuchos::OrdinalTraits<GO>::invalid() );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries-1), myImageID == numImages-1 ? 1 : Teuchos::OrdinalTraits<LO>::invalid() );
    TEST_COMPARE_ARRAYS( map.getNodeElementList(), myGlobal);
    TEST_EQUALITY_CONST( map.isNodeLocalElement(0), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(1), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(2), false ); // just try a couple
    TEST_EQUALITY_CONST( map.isNodeLocalElement(3), false );
    for (GO i=0; i < Teuchos::as<GO>(numGlobalEntries); ++i) {
      if (std::find(myGlobal.begin(),myGlobal.end(),i) == myGlobal.end()) {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), false );
      }
      else {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), true );
      }
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  //
  // INSTANTIATIONS
  //
#ifdef HAVE_XPETRA_TPETRA

  #define XPETRA_TPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N;

#endif

// List of tests (which run both on Epetra and Tpetra)
#define XP_MAP_INSTANT(LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, invalidConstructor1, M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, invalidConstructor2, M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, compatabilityTests,  M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, localMap,            M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, sameasTests,         M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, ContigUniformMap,    M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, validConstructor1,   M##LO##GO##N, LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, validConstructor2,   M##LO##GO##N, LO, GO, N)

// List of tests (which run on Tpetra only)
#define XPT_MAP_INSTANT(LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Map, validConstructor3,   M##LO##GO##N, LO, GO, N)


#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
TPETRA_INSTANTIATE_LGN ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_LGN ( XP_MAP_INSTANT )
TPETRA_INSTANTIATE_LGN ( XPT_MAP_INSTANT )

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(int,int,EpetraNode)
XP_MAP_INSTANT(int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(int,LongLong,EpetraNode)
XP_MAP_INSTANT(int,LongLong,EpetraNode)
#endif
#endif

}

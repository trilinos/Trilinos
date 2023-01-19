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

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_StridedMap.hpp"

#include "Xpetra_TpetraMap.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

namespace {

  bool testMpi = true;
  double errorTolSlack = 1e+1;

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
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, Constructor1, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 10 * numImages;
    std::vector<size_t> stridedInfo(1,1);
    SM map(lib, numGlobalElements, 0,stridedInfo, comm);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 1 );
    TEST_EQUALITY_CONST( map.isStrided(), false );
    TEST_EQUALITY_CONST( map.isBlocked(), false );

    stridedInfo.clear();
    stridedInfo.push_back(2);
    stridedInfo.push_back(1);
    SM map2(lib, 99, 0,stridedInfo, comm);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 3 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );

    stridedInfo.clear();
    stridedInfo.push_back(2);
    SM map3(lib, 100, 0,stridedInfo, comm);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 2 );
    TEST_EQUALITY_CONST( map3.isStrided(), false );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, Constructor2, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 10 * numImages;
    size_t numLocalElements = 10;
    std::vector<size_t> stridedInfo(1,1);

    SM map(lib, numGlobalElements, numLocalElements, 0, stridedInfo, comm);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 1 );
    TEST_EQUALITY_CONST( map.isStrided(), false );
    TEST_EQUALITY_CONST( map.isBlocked(), false );

    numGlobalElements = 33 * numImages;
    numLocalElements = 33;
    stridedInfo.clear();
    stridedInfo.push_back(2);
    stridedInfo.push_back(1);

    SM map2(lib, numGlobalElements, numLocalElements, 0, stridedInfo, comm);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 3 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );

    numGlobalElements = 20 * numImages;
    numLocalElements = 20;
    stridedInfo.clear();
    stridedInfo.push_back(2);
    SM map3(lib, numGlobalElements, numLocalElements, 0, stridedInfo, comm);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 2 );
    TEST_EQUALITY_CONST( map3.isStrided(), false );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, Constructor3, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // constructor calls: (num global elements, index base)
    GO indexBase = 1111;
    Xpetra::global_size_t numGlobalElements = 10 * numImages;
    size_t numLocalElements = 10;
    std::vector<size_t> stridedInfo(1,1);

    SM map(lib, numGlobalElements, numLocalElements, indexBase, stridedInfo, comm);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 1 );
    TEST_EQUALITY_CONST( map.isStrided(), false );
    TEST_EQUALITY_CONST( map.isBlocked(), false );
    TEST_EQUALITY_CONST( map.getIndexBase(), indexBase );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), indexBase + Teuchos::as<GO>(numGlobalElements) - 1);

    numGlobalElements = 33 * numImages;
    numLocalElements = 33;
    stridedInfo.clear();
    stridedInfo.push_back(2);
    stridedInfo.push_back(1);

    SM map2(lib, numGlobalElements, numLocalElements, indexBase, stridedInfo, comm);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 3 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getIndexBase(), indexBase );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), indexBase + Teuchos::as<GO>(numGlobalElements) - 1);

    numGlobalElements = 20 * numImages;
    numLocalElements = 20;
    stridedInfo.clear();
    stridedInfo.push_back(2);
    SM map3(lib, numGlobalElements, numLocalElements, indexBase, stridedInfo, comm);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 2 );
    TEST_EQUALITY_CONST( map3.isStrided(), false );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
    TEST_EQUALITY_CONST( map3.getIndexBase(), indexBase );
    TEST_EQUALITY_CONST( map3.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST( map3.getMaxAllGlobalIndex(), indexBase + Teuchos::as<GO>(numGlobalElements) - 1);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, StridedPartConstructor1, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 30 * numImages;
    size_t numLocalElements = 30;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(2);
    stridedInfo.push_back(1);

    SM map(lib, numGlobalElements, numLocalElements, 0,stridedInfo, comm, 0);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 3 );
    TEST_EQUALITY_CONST( map.isStrided(), true );
    TEST_EQUALITY_CONST( map.isBlocked(), true );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), 0 );
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) - 2 );
    TEST_EQUALITY_CONST( map.isContiguous(), false);
    TEST_EQUALITY_CONST( map.getLocalNumElements() % 2 , 0);
    TEST_EQUALITY_CONST( map.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[0] );

    SM map2(lib, numGlobalElements, numLocalElements, 0,stridedInfo, comm, 1);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 3 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), 2 );
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map2.isContiguous(), false);
    TEST_EQUALITY_CONST( map2.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[1]);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, StridedPartConstructor2, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    SM map(lib, numGlobalElements, 0,stridedInfo, comm, 0);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map.isStrided(), true );
    TEST_EQUALITY_CONST( map.isBlocked(), true );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), 0 );
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) - 10 );
    TEST_EQUALITY_CONST( map.isContiguous(), false);
    TEST_EQUALITY_CONST( map.getLocalNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[0] );

    SM map2(lib, numGlobalElements, 0,stridedInfo, comm, 1);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), 3 );
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) - 6 );
    TEST_EQUALITY_CONST( map2.isContiguous(), false);
    TEST_EQUALITY_CONST( map2.getLocalNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map2.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[1]);

    SM map3(lib, numGlobalElements, 0,stridedInfo, comm, 2);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3.isStrided(), true );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
    TEST_EQUALITY_CONST( map3.getMinAllGlobalIndex(), 7 );
    TEST_EQUALITY_CONST( map3.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map3.isContiguous(), false);
    TEST_EQUALITY_CONST( map3.getLocalNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map3.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[2]);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, StridedPartConstructor3, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    // indexBase = 1111
    GO indexBase = 1111;

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    SM map(lib, numGlobalElements, indexBase ,stridedInfo, comm, 0);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map.isStrided(), true );
    TEST_EQUALITY_CONST( map.isBlocked(), true );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), indexBase );
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) + indexBase - 10 );
    TEST_EQUALITY_CONST( map.isContiguous(), false);
    TEST_EQUALITY_CONST( map.getLocalNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[0] );

    SM map2(lib, numGlobalElements, indexBase ,stridedInfo, comm, 1);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), 3 + indexBase );
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) + indexBase - 6 );
    TEST_EQUALITY_CONST( map2.isContiguous(), false);
    TEST_EQUALITY_CONST( map2.getLocalNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map2.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[1]);

    SM map3(lib, numGlobalElements, indexBase ,stridedInfo, comm, 2);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3.isStrided(), true );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
    TEST_EQUALITY_CONST( map3.getMinAllGlobalIndex(), 7 + indexBase);
    TEST_EQUALITY_CONST( map3.getMaxAllGlobalIndex(), Teuchos::as<GO>(numGlobalElements) + indexBase - 1 );
    TEST_EQUALITY_CONST( map3.isContiguous(), false);
    TEST_EQUALITY_CONST( map3.getLocalNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map3.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[2]);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, StridedPartConstructorWithOffset, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    GO offset = 111;

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    SM map(lib, numGlobalElements, 0,stridedInfo, comm, 0, offset);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map.isStrided(), true );
    TEST_EQUALITY_CONST( map.isBlocked(), true );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), offset );
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 10 );
    TEST_EQUALITY_CONST( map.isContiguous(), false);
    TEST_EQUALITY_CONST( map.getLocalNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[0] );

    SM map2(lib, numGlobalElements, 0,stridedInfo, comm, 1, offset);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), offset + 3 );
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 6 );
    TEST_EQUALITY_CONST( map2.isContiguous(), false);
    TEST_EQUALITY_CONST( map2.getLocalNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map2.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[1]);

    SM map3(lib, numGlobalElements, 0,stridedInfo, comm, 2, offset);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3.isStrided(), true );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
    TEST_EQUALITY_CONST( map3.getMinAllGlobalIndex(), offset + 7 );
    TEST_EQUALITY_CONST( map3.getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map3.isContiguous(), false);
    TEST_EQUALITY_CONST( map3.getLocalNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map3.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[2]);
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( StridedMap, StridedPartConstructorOffsetPlusIndexBase, M, LO, GO, N )
  {
    typedef Xpetra::StridedMap<LO,GO,N> SM;

    // create a comm
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // test constructor for Xpetra::StridedMaps
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    GO offset = 111;
    GO indexBase = 89;

    // constructor calls: (num global elements, index base)
    Xpetra::global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    SM map(lib, numGlobalElements, indexBase, stridedInfo, comm, 0, offset);
    TEST_EQUALITY_CONST( map.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map.isStrided(), true );
    TEST_EQUALITY_CONST( map.isBlocked(), true );
    TEST_EQUALITY_CONST( map.getMinAllGlobalIndex(), indexBase + offset );
    TEST_EQUALITY_CONST( map.getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 10 );
    TEST_EQUALITY_CONST( map.isContiguous(), false);
    TEST_EQUALITY_CONST( map.getLocalNumElements() % 3 , 0);
    TEST_EQUALITY_CONST( map.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[0] );

    SM map2(lib, numGlobalElements, indexBase ,stridedInfo, comm, 1, offset);
    TEST_EQUALITY_CONST( map2.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map2.isStrided(), true );
    TEST_EQUALITY_CONST( map2.isBlocked(), true );
    TEST_EQUALITY_CONST( map2.getMinAllGlobalIndex(), indexBase + offset + 3 );
    TEST_EQUALITY_CONST( map2.getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 6 );
    TEST_EQUALITY_CONST( map2.isContiguous(), false);
    TEST_EQUALITY_CONST( map2.getLocalNumElements() % 4 , 0);
    TEST_EQUALITY_CONST( map2.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[1]);

    SM map3(lib, numGlobalElements, indexBase ,stridedInfo, comm, 2, offset);
    TEST_EQUALITY_CONST( map3.getFixedBlockSize(), 12 );
    TEST_EQUALITY_CONST( map3.isStrided(), true );
    TEST_EQUALITY_CONST( map3.isBlocked(), true );
    TEST_EQUALITY_CONST( map3.getMinAllGlobalIndex(), indexBase + offset + 7 );
    TEST_EQUALITY_CONST( map3.getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 1 );
    TEST_EQUALITY_CONST( map3.isContiguous(), false);
    TEST_EQUALITY_CONST( map3.getLocalNumElements() % 5 , 0);
    TEST_EQUALITY_CONST( map3.getLocalNumElements(), numLocalElements / map.getFixedBlockSize() * stridedInfo[2]);
  }


  //
  // INSTANTIATIONS
  //

  #define XPETRA_TPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N;

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N;

#endif

// List of tests (which run both on Epetra and Tpetra)
#define XP_MAP_INSTANT(LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, Constructor1, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, Constructor2, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, Constructor3, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, StridedPartConstructor1, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, StridedPartConstructor2, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, StridedPartConstructor3, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, StridedPartConstructorWithOffset, M##LO##GO##N , LO, GO, N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( StridedMap, StridedPartConstructorOffsetPlusIndexBase, M##LO##GO##N , LO, GO, N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
TPETRA_INSTANTIATE_LGN ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_LGN ( XP_MAP_INSTANT )

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

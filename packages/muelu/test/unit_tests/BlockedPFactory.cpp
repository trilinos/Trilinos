// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  /////////////////////////
  // helper function

  Teuchos::RCP<CrsMatrixWrap> GenerateProblemMatrix(const Teuchos::RCP<const Map> rangemap, const Teuchos::RCP<const Map> domainmap, Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {

    Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map,CrsMatrixWrap>::Build(rangemap, 3);

    LocalOrdinal NumMyRowElements = rangemap->getNodeNumElements();

    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalColElements = domainmap->getNodeElementList();
    GlobalOrdinal NumGlobalColElements = domainmap->getGlobalNumElements();
    //GlobalOrdinal nIndexBase = domainmap->getIndexBase();

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz=2;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    for (LocalOrdinal i = 0; i < NumMyRowElements; ++i)
    {
      if(i < MyGlobalColElements.size()) {
        if (MyGlobalColElements[i] == domainmap->getMinGlobalIndex())
        {
          // off-diagonal for first row
          Indices[0] = domainmap->getMinGlobalIndex();
          NumEntries = 1;
          Values[0] = c;
        }
        else if (MyGlobalColElements[i] == domainmap->getMinGlobalIndex() + NumGlobalColElements - 1)
        {
          // off-diagonal for last row
          Indices[0] = domainmap->getMinGlobalIndex() + NumGlobalColElements - 2;
          NumEntries = 1;
          Values[0] = b;
        }
        else
        {
          // off-diagonal for internal row
          Indices[0] = MyGlobalColElements[i] - 1;
          Values[1] = b;
          Indices[1] = MyGlobalColElements[i] + 1;
          Values[0] = c;
          NumEntries = 2;
        }

        // put the off-diagonal entries
        // Xpetra wants ArrayViews (sigh)
        Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
        Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
        mtx->insertGlobalValues(rangemap->getGlobalElement(i), iv, av);

        // Put in the diagonal entry
        mtx->insertGlobalValues(rangemap->getGlobalElement(i),
            Teuchos::tuple<GlobalOrdinal>(MyGlobalColElements[i]),
            Teuchos::tuple<Scalar>(a) );

      }

    } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)

    mtx->fillComplete(domainmap,rangemap);

    return mtx;
  }

  TEUCHOS_UNIT_TEST(BlockedPFactory, Constructor)
  {
    // test for accessing subblocks from a blocked CRS Matrix using SubBlockAFactory
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    RCP<const Map> bigMap;
    RCP<const Map> map1;
    RCP<const Map> map2;
    GO numElements = 400;
    GO numElements1 = 200;
    GO numElements2 = 200;

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(1);

    map1   = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);
    map2   = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm);

    std::vector<GlobalOrdinal> localGids; // vector with all local GIDs on cur proc
    Teuchos::ArrayView< const GlobalOrdinal > map1eleList = map1->getNodeElementList(); // append all local gids from map1 and map2
    localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
    Teuchos::ArrayView< const GlobalOrdinal > map2eleList = map2->getNodeElementList();
    localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
    Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0],localGids.size());
    bigMap = MapFactory::Build(lib, numElements, eleList, 0, comm); // create full big map (concatenation of map1 and map2)

    std::vector<Teuchos::RCP<const Map> > maps;
    maps.push_back(map1); maps.push_back(map2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>::Build(bigMap, maps);

    RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix(map1,map1,2,-1,-1);
    RCP<CrsMatrixWrap> Op12 = GenerateProblemMatrix(map1,map2,1, 0, 0);
    RCP<CrsMatrixWrap> Op21 = GenerateProblemMatrix(map2,map1,1, 0, 0);
    RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix(map2,map2,3,-2,-1);

    //Op11->describe(out, Teuchos::VERB_EXTREME);

    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node,LocalMatOps>(mapExtractor,mapExtractor,10));

    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat11 = Op11->getCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat12 = Op12->getCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat21 = Op21->getCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat22 = Op22->getCrsMatrix();
    bOp->setMatrix(0,0,crsMat11);
    bOp->setMatrix(0,1,crsMat12);
    bOp->setMatrix(1,0,crsMat21);
    bOp->setMatrix(1,1,crsMat22);
    bOp->fillComplete();
    TEST_EQUALITY(bOp!=Teuchos::null, true);

    // build hierarchy
    RCP<Level> levelOne = rcp(new Level());
    RCP<Level> levelTwo = rcp(new Level()); levelTwo->SetPreviousLevel(levelOne);
    levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp)); // set blocked operator

    // define sub block factories for blocked operator "A"
    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    RCP<SubBlockAFactory> A12Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 1));
    RCP<SubBlockAFactory> A21Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 0));
    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    // request subblocks of A
    levelOne->Request("A", A11Fact.get(), MueLu::NoFactory::get());
    levelOne->Request("A", A22Fact.get(), MueLu::NoFactory::get());
    TEST_EQUALITY(levelOne->IsRequested("A", A11Fact.get()),true);
    TEST_EQUALITY(levelOne->IsRequested("A", A22Fact.get()),true);

    RCP<Matrix> A11 = levelOne->Get<RCP<Matrix> >("A",A11Fact.get());
    RCP<Matrix> A22 = levelOne->Get<RCP<Matrix> >("A",A22Fact.get());
    TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()),true);
    TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()),true);

    // store data A11 and A22 as variable "P" on level 2

    levelTwo->Request("P", A11Fact.get(), MueLu::NoFactory::get());
    levelTwo->Request("P", A22Fact.get(), MueLu::NoFactory::get());
    TEST_EQUALITY(levelTwo->IsRequested("P", A11Fact.get()),true);
    TEST_EQUALITY(levelTwo->IsRequested("P", A22Fact.get()),true);

    levelTwo-> Set("P", A11, A11Fact.get());
    levelTwo-> Set("P", A22, A22Fact.get());
    TEST_EQUALITY(levelTwo->IsAvailable("P", A11Fact.get()),true);
    TEST_EQUALITY(levelTwo->IsAvailable("P", A22Fact.get()),true);

    levelOne->Release("A", A11Fact.get());
    levelOne->Release("A", A22Fact.get());

    // set up factory manager
    RCP<FactoryManager> FC1 = rcp(new FactoryManager());
    FC1->SetFactory("P", A11Fact);  // fool P to be generated by A11Fact
    FC1->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    RCP<FactoryManager> FC2 = rcp(new FactoryManager());
    FC2->SetFactory("P", A22Fact);  // fool P to be generated by A11Fact
    FC2->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

    /////////////////////////////////////////// define blocked transfer ops
    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory(Teuchos::null)); // use row map index base from bOp
    PFact->AddFactoryManager(FC1);
    PFact->AddFactoryManager(FC2);

    levelTwo->Request("P", PFact.get(), MueLu::NoFactory::get());
    TEST_EQUALITY(levelTwo->IsRequested("P", PFact.get()),true);

    RCP<Matrix> P = levelTwo->Get<RCP<Matrix> >("P",PFact.get());
    TEST_EQUALITY(P!=Teuchos::null,true);

    RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > bP = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node,LocalMatOps> >(P);
    TEST_EQUALITY(bP!=Teuchos::null,true);

    TEST_EQUALITY(bOp->Rows(),2);
    TEST_EQUALITY(bOp->Cols(),2);

    // TODO add some more tests
  } //Constructor
}



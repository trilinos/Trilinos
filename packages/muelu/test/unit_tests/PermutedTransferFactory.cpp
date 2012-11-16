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
#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_PermutedTransferFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#include "Galeri_XpetraUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(PermutedTransfer, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<PermutedTransferFactory> ptFactory = rcp(new PermutedTransferFactory);
    TEST_EQUALITY(ptFactory != Teuchos::null, true);
  } // Constructor test

  TEUCHOS_UNIT_TEST(PermutedTransfer, Build1)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    GO nx = 199;
    RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx);
    fineLevel.Set("A",A);

    //build coordinates
    Teuchos::ParameterList list;
    list.set("nx",nx);
    RCP<MultiVector> coordVector = Galeri::Xpetra::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",A->getRowMap(),list);
    fineLevel.Set("Coordinates",coordVector);


    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    RCP<TentativePFactory>    Ptentfact = rcp(new TentativePFactory(UCAggFact));
    RCP<SaPFactory>           Pfact = rcp( new SaPFactory(Ptentfact));
    RCP<FactoryBase2>             Rfact = rcp( new TransPFactory(Pfact) );
    RCP<RAPFactory>           Acfact = rcp( new RAPFactory(Pfact,Rfact) );
    RCP<FactoryBase2>             Rtentfact = rcp( new TransPFactory(Ptentfact) );

    RCP<MultiVectorTransferFactory> mvTransFact = rcp(new MultiVectorTransferFactory("Coordinates","R",Rtentfact));
    Acfact->AddTransferFactory(mvTransFact);
    RCP<ZoltanInterface>      zoltan = rcp(new ZoltanInterface(Acfact,mvTransFact));
    RCP<RepartitionFactory> RepartitionFact = rcp(new RepartitionFactory(zoltan,Acfact));

    coarseLevel.Request("A",Acfact.get());  // kick off the DeclareInputs
    coarseLevel.Request("Permutation",RepartitionFact.get());  // request permutation matrix
    //coarseLevel.Request("P",Pfact.get());
    coarseLevel.Request("R",Rtentfact.get());
    coarseLevel.Request("Coordinates",mvTransFact.get());

    RCP<PermutedTransferFactory> ptFactory = rcp( new PermutedTransferFactory(RepartitionFact, Acfact, Pfact, MueLu::INTERPOLATION) );
    coarseLevel.Request("P",ptFactory.get());
    ptFactory->Build(fineLevel,coarseLevel);

    ptFactory = rcp( new PermutedTransferFactory(RepartitionFact, Acfact, Rfact, MueLu::RESTRICTION) );
    coarseLevel.Request("R",ptFactory.get());
    ptFactory->Build(fineLevel,coarseLevel);



  } // Constructor test

} // namespace MueLuTests


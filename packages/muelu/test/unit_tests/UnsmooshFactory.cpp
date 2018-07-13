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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_VariableDofLaplacianFactory.hpp"
#include "MueLu_UnsmooshFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_HierarchyUtils.hpp"

namespace MueLuTests {


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UnsmooshFactory, UnsmooshTentativeP, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    if (!TYPE_EQUAL(GO, int)) { out << "Skipping test for GO != int"        << std::endl; return; }
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // TAW 04/21: test is crashing on 4 procs with Epetra due to an unknown reason in the Epetra_BlockMap constructor (MPI communication)
    if ( comm->getSize() > 2 && lib == Xpetra::UseEpetra) { out << "Skipping test for more than 2 procs when using Epetra"        << std::endl; return; }

    GlobalOrdinal nx = 6, ny = 6;

    typedef Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node> mv_type_double;
    typedef Xpetra::MultiVectorFactory<double,LocalOrdinal,GlobalOrdinal,Node> MVFactory_double;

    // Describes the initial layout of matrix rows across processors.
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(lib, "Cartesian2D", comm, galeriList);

    //build coordinates before expanding map (nodal coordinates, not dof-based)
    RCP<mv_type_double> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LocalOrdinal,GlobalOrdinal,Map,mv_type_double>("2D", nodeMap, galeriList);
    RCP<const Map> dofMap = MapFactory::Build(nodeMap, 2); //expand map for 2 DOFs per node

    galeriList.set("right boundary" , "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary"   , "Neumann");
    galeriList.set("front boundary" , "Neumann");
    galeriList.set("back boundary"  , "Neumann");
    galeriList.set("keepBCs",             false);

    int maxDofPerNode = 2;

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", dofMap, galeriList);
    RCP<Matrix> A = Pr->BuildMatrix();
    A->SetFixedBlockSize(maxDofPerNode);

    //->describe(out, Teuchos::VERB_EXTREME);

    TEST_EQUALITY(dofMap->getNodeNumElements(),2*nodeMap->getNodeNumElements());

    // TODO introduce extra factory manager objects


    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    // Test of createTwoLevelHierarchy: to be moved...
    TEST_EQUALITY(fineLevel.GetLevelID(), 0);
    TEST_EQUALITY(coarseLevel.GetLevelID(), 1);



    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates",coordinates);

    // DofPresent for VariableDofLaplacianFactory
    Teuchos::ArrayRCP<LocalOrdinal> dofPresent(A->getRowMap()->getNodeNumElements(),1);
    fineLevel.Set<Teuchos::ArrayRCP<LocalOrdinal> >("DofPresent", dofPresent);
    // DofStatus for UnsmooshFactory
    Teuchos::Array<char> dofStatus = Teuchos::Array<char>(A->getRangeMap()->getNodeNumElements() * maxDofPerNode,'s');
    fineLevel.Set("DofStatus", dofStatus);

    VariableDofLaplacianFactory lapFact;
    lapFact.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(maxDofPerNode));

    AmalgamationFactory amalgFact;
    CoalesceDropFactory coalFact;
    NullspaceFactory nspFact;
    TentativePFactory Pfact;
    UnsmooshFactory unsmooFact;
    FactoryManager innerFactManager;
    innerFactManager.SetKokkosRefactor(false);
    innerFactManager.SetFactory("A", Teuchos::rcpFromRef(lapFact));
    innerFactManager.SetFactory("UnAmalgamationInfo", Teuchos::rcpFromRef(amalgFact));
    innerFactManager.SetFactory("Graph", Teuchos::rcpFromRef(coalFact));
    innerFactManager.SetFactory("DofsPerNode", Teuchos::rcpFromRef(coalFact));
    innerFactManager.SetFactory("Nullspace", Teuchos::rcpFromRef(nspFact));
    //innerFactManager.SetIgnoreUserData(true);

    fineLevel.SetFactoryManager(Teuchos::rcpFromRef(innerFactManager));
    coarseLevel.SetFactoryManager(Teuchos::rcpFromRef(innerFactManager));

    // TODO it would be nice if we could use an innerFactoryManager for building Ptent and use
    //      it for the outer factory manager.
    //{
    //  MueLu::SetFactoryManager fineSFM  (Teuchos::rcpFromRef(fineLevel),   Teuchos::rcpFromRef(innerFactManager));
    //  MueLu::SetFactoryManager coarseSFM(Teuchos::rcpFromRef(coarseLevel), Teuchos::rcpFromRef(innerFactManager));
      nspFact.SetFactory("A", Teuchos::rcpFromRef(lapFact));   // we have to provide the correct A (since we cannot set ignoreUserData
      amalgFact.SetFactory("A", Teuchos::rcpFromRef(lapFact));
      coalFact.SetFactory("A", Teuchos::rcpFromRef(lapFact));
      Pfact.SetFactory("A", Teuchos::rcpFromRef(lapFact));
    //}

    unsmooFact.SetFactory("P", Teuchos::rcpFromRef(Pfact));
    unsmooFact.SetFactory("DofStatus", MueLu::NoFactory::getRCP());
    unsmooFact.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(maxDofPerNode));
    unsmooFact.SetParameter("fineIsPadded", Teuchos::ParameterEntry(true));
    unsmooFact.SetFactory("A", MueLu::NoFactory::getRCP());
    coarseLevel.Request("P",&unsmooFact);
    RCP<Matrix> Pfinal = coarseLevel.Get<RCP<Matrix> >("P",&unsmooFact);

    //coarseLevel.print(out,MueLu::Extreme);
    //Pfinal->describe(out, Teuchos::VERB_EXTREME);

    TEST_EQUALITY(Pfinal->getRowMap()->isSameAs(*A->getRowMap()),true);
  } // UnsmooshTentativeP

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UnsmooshFactory, UnsmooshTentativeP, SC, LO, GO, Node) \


#include <MueLu_ETI_4arg.hpp>
}

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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <stdlib.h>

#include <Teuchos_UnitTestHarness.hpp>
// #include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_NoFactory.hpp"

#include "MueLu_CreateXpetraPreconditioner.hpp"

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_SERIAL)

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UserData, CreateXpetraPreconditioner, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using Teuchos::RCP;

  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  GO nx = 1000;

  std::string xmlFileName                          = "UserData/test.xml";
  Teuchos::RCP<Teuchos::ParameterList> inParamList = Teuchos::getParametersFromXmlFile(xmlFileName);
  if (lib == Xpetra::UseEpetra) {
    inParamList->sublist("Hierarchy").set("use kokkos refactor", false);
  }

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

  // Matrix
  RCP<Matrix> Op     = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx * comm->getSize(), lib);
  RCP<const Map> map = Op->getRowMap();

  // Normalized RHS
  RCP<MultiVector> RHS1 = MultiVectorFactory::Build(map, 1);
  RHS1->setSeed(846930886);
  RHS1->randomize();
  Teuchos::Array<magnitude_type> norms(1);
  RHS1->norm2(norms);
  RHS1->scale(1 / norms[0]);

  // Zero initial guess
  RCP<MultiVector> X1 = MultiVectorFactory::Build(Op->getRowMap(), 1);
  X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", Op->getRowMap(), galeriList);
  RCP<MultiVector> nullspace   = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(Op->getDomainMap(), 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

  // Add sublist "user data" to MueLu's parameter list
  const std::string userName            = "user data";
  Teuchos::ParameterList& userParamList = inParamList->sublist(userName);

  // Create test variables to be stored on Level 0 of the Hierarchy
  SC myScalar = 3.14;
  userParamList.set<SC>("Scalar myScalar", myScalar);
  double myDouble = 2.71;
  userParamList.set<double>("Double myDouble", myDouble);
  int myInt = 42;
  userParamList.set<int>("int myInt", myInt);
  std::string myString = "Test string";
  userParamList.set<std::string>("String myString", myString);
  Array<GO> myArrayGO(4);
  myArrayGO[0] = 1;
  myArrayGO[1] = 4;
  myArrayGO[2] = 5;
  myArrayGO[3] = 0;
  userParamList.set<Array<GO> >("Array<GO> myArray<GO>", myArrayGO);
  Array<LO> myArrayLO(5);
  myArrayLO[0] = 8;
  myArrayLO[1] = 7;
  myArrayLO[2] = 1;
  myArrayLO[3] = 2;
  myArrayLO[4] = 3;
  userParamList.set<Array<LO> >("Array<LO> myArray<LO>", myArrayLO);

  RCP<Hierarchy> xH = MueLu::CreateXpetraPreconditioner<SC, LO, GO, NO>(Op, *inParamList);

  // Extract variables on level 0 and check that they are unchanged.
  RCP<MueLu::Level> level0 = xH->GetLevel();
  bool result              = true;
  std::string errorMsg     = "";

  if (!(level0->Get<SC>("myScalar") == myScalar)) {
    errorMsg += "myScalar does not have correct value on level 0.\n";
    result = false;
  }

  if (!(level0->Get<double>("myDouble") == myDouble)) {
    errorMsg += "myDouble does not have correct value on level 0.\n";
    result = false;
  }

  if (!(level0->Get<int>("myInt") == myInt)) {
    errorMsg += "myInt does not have correct value on level 0.\n";
    result = false;
  }

  if (!(level0->Get<std::string>("myString") == myString)) {
    errorMsg += "myString does not have correct value on level 0.\n";
    result = false;
  }

  if (!(level0->Get<Array<GO> >("myArray<GO>") == myArrayGO)) {
    errorMsg += "myArray<GO> does not have correct value on level 0.\n";
    result = false;
  }

  if (!(level0->Get<Array<LO> >("myArray<LO>") == myArrayLO)) {
    errorMsg += "myArray<LO> does not have correct value on level 0.\n";
    result = false;
  }

  TEST_EQUALITY(result, true);
}  // CreatePreconditioner

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UserData, CreateXpetraPreconditioner, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

#endif  // HAVE_MUELU_EPETRA && HAVE_MUELU_SERIAL

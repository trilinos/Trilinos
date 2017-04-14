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

#include "MueLu_VariableDofLaplacianFactory.hpp"

namespace MueLuTests {


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    GlobalOrdinal nEle = 63;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",nEle);
    matrixParameters.set("ny",nEle);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,CrsMatrixWrap,MultiVector>("Laplace2D", map, matrixParameters);
    RCP<Matrix> A = Pr->BuildMatrix();
    RCP<MultiVector> coords = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", A->getRowMap(), matrixParameters);

    // build hierarchy
    RCP<Level> l = rcp(new Level());
    l->SetLevelID(0);
    l->SetComm(comm);
    l->Set("A", A);
    l->Set("Coordinates",coords);

    Teuchos::ArrayRCP<const bool> dofPresent(A->getRowMap()->getNodeNumElements(),true);
    l->Set<Teuchos::ArrayRCP<const bool> >("DofPresent", dofPresent);

    VariableDofLaplacianFactory lapFact;

    l->Request("A",&lapFact);

    l->Get<RCP<Matrix> >("A",&lapFact);


    //TEST_EQUALITY(AA.get(), A.get());
  } // InputData


  // helper class for testing functionality of FineLevelInputDataFactory
  /*template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class FineLevelInputDataFactoryTester {
#   include "MueLu_UseShortNames.hpp"
  public:
    void test_testfunc(const FineLevelInputDataFactory& fac) {
      std::cout << "FineLevelInputDataFactoryTester" << std::endl;
      fac.test();
    }
  };*/

  // unit test just for demonstration purposes
  /*TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FineLevelInputDataFactory, TestFunc, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    FineLevelInputDataFactory fac;

    FineLevelInputDataFactoryTester<Scalar,LocalOrdinal,GlobalOrdinal,Node> tester;

    tester.test_testfunc(fac);

    //TEST_EQUALITY(AA.get(), A.get());
  } // TestFunc*/

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(VariableDofLaplacianFactory, VarLaplConstructor, SC, LO, GO, Node) \

#include <MueLu_ETI_4arg.hpp>
}



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

// declare content from Galeri_XpetraMaps.hpp
// we cannot include the header file, since it
// is already included for the Repartition.cpp
// unit tests
//#include <Galeri_XpetraMaps.hpp>
namespace Galeri {
  namespace Xpetra {

    using Teuchos::RCP;

    //! Map creation function (for Tpetra, Epetra, Xpetra::TpetraMap and Xpetra::EpetraMap)
    template <class LocalOrdinal, class GlobalOrdinal, class Map>
    RCP<Map> CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList & list);

#ifdef HAVE_GALERI_XPETRA
    //! Map creation function (for Xpetra::Map with an UnderlyingLib parameter)
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList & list);
#endif
  }
}

#include "MueLu_VariableDofLaplacianFactory.hpp"

namespace MueLuTests {


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    if (!TYPE_EQUAL(GO, int)) { out << "Skipping test for GO != int"        << std::endl; return; }
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    GlobalOrdinal nEle = 64;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
    Teuchos::ParameterList matrixParameters;
    matrixParameters.set("nx",8);
    matrixParameters.set("ny",8);

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

    RCP<Matrix> lapA = l->Get<RCP<Matrix> >("A",&lapFact);

    //lapA->describe(out, Teuchos::VERB_EXTREME);

    TEST_EQUALITY(lapA->getRowMap()->isSameAs(*A->getRowMap()),true);

    Teuchos::RCP<Vector> oneVec = VectorFactory::Build(A->getRowMap());
    Teuchos::RCP<Vector> res = VectorFactory::Build(A->getRowMap());
    oneVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    res->putScalar(27*Teuchos::ScalarTraits<Scalar>::one());
    Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(lapA)->apply(*oneVec,*res);
    TEST_COMPARE(res->normInf(),<, 1e-13);
  } // VarLaplConstructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplConstructor2, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    if (!TYPE_EQUAL(GO, int)) { out << "Skipping test for GO != int"        << std::endl; return; }
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    GlobalOrdinal nx = 4, ny = 4;

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

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", dofMap, galeriList);
    RCP<Matrix> A = Pr->BuildMatrix();
    A->SetFixedBlockSize(2);

    //->describe(out, Teuchos::VERB_EXTREME);

    TEST_EQUALITY(dofMap->getNodeNumElements(),2*nodeMap->getNodeNumElements());

    // build hierarchy
    RCP<Level> l = rcp(new Level());
    l->SetLevelID(0);
    l->SetComm(comm);
    l->Set("A", A);
    l->Set("Coordinates",coordinates);

    Teuchos::ArrayRCP<const bool> dofPresent(A->getRowMap()->getNodeNumElements(),true);
    l->Set<Teuchos::ArrayRCP<const bool> >("DofPresent", dofPresent);

    A->getColMap()->describe(out,Teuchos::VERB_EXTREME);

    VariableDofLaplacianFactory lapFact;
    lapFact.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(2));
    l->Request("A",&lapFact);

    RCP<Matrix> lapA = l->Get<RCP<Matrix> >("A",&lapFact);

    lapA->describe(out, Teuchos::VERB_EXTREME);

    /*Teuchos::RCP<Vector> oneVec = VectorFactory::Build(lapA->getRowMap());
    Teuchos::RCP<Vector> res = VectorFactory::Build(lapA->getRowMap());
    oneVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    res->putScalar(27*Teuchos::ScalarTraits<Scalar>::one());
    Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(lapA)->apply(*oneVec,*res);
    TEST_COMPARE(res->normInf(),<, 1e-13);*/
  } // VarLaplConstructor2

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
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(VariableDofLaplacianFactory, VarLaplConstructor2, SC, LO, GO, Node) \

#include <MueLu_ETI_4arg.hpp>
}



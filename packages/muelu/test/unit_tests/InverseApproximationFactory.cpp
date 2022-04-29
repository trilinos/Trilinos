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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include <MueLu_InverseApproximationFactory.hpp>

namespace MueLuTests {

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InverseApproximationFactory, InverseDiagonalConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    {
      const int n = 20;
      Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build1DPoisson(n, lib);

      Level level;
      TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<InverseApproximationFactory> invapproxFact = rcp( new InverseApproximationFactory() );
      invapproxFact->SetFactory("A",MueLu::NoFactory::getRCP());
      invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("diagonal")));

      // request InverseApproximation operator
      level.Request("Ainv", invapproxFact.get());

      // generate Schur complement operator
      invapproxFact->Build(level);

      RCP<Vector> Ainv = level.Get<RCP<Vector> >("Ainv", invapproxFact.get());
      TEST_EQUALITY(Ainv.is_null(), false);
      TEST_EQUALITY(Ainv->getMap()->getMinGlobalIndex(), comm->getRank() * int(n/comm->getSize()));
      TEST_EQUALITY(Ainv->getMap()->getMaxGlobalIndex(), comm->getRank() * int(n/comm->getSize()) + int(n/comm->getSize()-1));

      Teuchos::ArrayRCP<const Scalar> AinvData = Ainv->getData(0);
      bool bCheck = true;
      for(int i=0; i<int(n/comm->getSize()); i++) if(AinvData[i] != Teuchos::as<Scalar>(0.5)) bCheck = false;
      TEST_EQUALITY(bCheck, true);
    }

    {
      const int n = 42;
      Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build2DPoisson(n, n, lib);

      Level level;
      TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<InverseApproximationFactory> invapproxFact = rcp( new InverseApproximationFactory() );
      invapproxFact->SetFactory("A",MueLu::NoFactory::getRCP());
      invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("diagonal")));

      // request InverseApproximation operator
      level.Request("Ainv", invapproxFact.get());

      // generate Schur complement operator
      invapproxFact->Build(level);

      RCP<Vector> Ainv = level.Get<RCP<Vector> >("Ainv", invapproxFact.get());
      TEST_EQUALITY(Ainv.is_null(), false);
      TEST_EQUALITY(Ainv->getMap()->getMinGlobalIndex(), comm->getRank() * int(n*n/comm->getSize()));
      TEST_EQUALITY(Ainv->getMap()->getMaxGlobalIndex(), comm->getRank() * int(n*n/comm->getSize()) + int(n*n/comm->getSize()-1));

      Teuchos::ArrayRCP<const Scalar> AinvData = Ainv->getData(0);
      bool bCheck = true;
      for(int i=0; i<int(n*n/comm->getSize()); i++) if(AinvData[i] != Teuchos::as<Scalar>(0.25)) bCheck = false;
      TEST_EQUALITY(bCheck, true);
    }

  } //InverseDiagonalConstructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InverseApproximationFactory, InverseLumpingConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    {
      const int n = 20;
      Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build1DPoisson(n, lib);
      A->scale(0.5);

      Level level;
      TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<InverseApproximationFactory> invapproxFact = rcp( new InverseApproximationFactory() );
      invapproxFact->SetFactory("A",MueLu::NoFactory::getRCP());
      invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("lumping")));

      // request InverseApproximation operator
      level.Request("Ainv", invapproxFact.get());

      // generate Schur complement operator
      invapproxFact->Build(level);

      RCP<Vector> Ainv = level.Get<RCP<Vector> >("Ainv", invapproxFact.get());
      TEST_EQUALITY(Ainv.is_null(), false);
      TEST_EQUALITY(Ainv->getMap()->getMinGlobalIndex(), comm->getRank() * int(n/comm->getSize()));
      TEST_EQUALITY(Ainv->getMap()->getMaxGlobalIndex(), comm->getRank() * int(n/comm->getSize()) + int(n/comm->getSize()-1));

      Teuchos::ArrayRCP<const Scalar> AinvData = Ainv->getData(0);
      bool bCheck = false;
      for(int i=0; i<int(n/comm->getSize()); i++)
        if(std::abs(AinvData[i] - Teuchos::as<Scalar>(0.66666666667)) < std::abs(Teuchos::as<Scalar>(1e-8)) ||
           std::abs(AinvData[i] - Teuchos::as<Scalar>(0.5)) < std::abs(Teuchos::as<Scalar>(1e-8)))
          bCheck = true;
      TEST_EQUALITY(bCheck, true);
    }

    {
      const int n = 42;
      Teuchos::RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build2DPoisson(n, n, lib);

      Level level;
      TestHelpers::TestFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createSingleLevelHierarchy(level);
      level.Set("A", A);

      RCP<InverseApproximationFactory> invapproxFact = rcp( new InverseApproximationFactory() );
      invapproxFact->SetFactory("A",MueLu::NoFactory::getRCP());
      invapproxFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("lumping")));

      // request InverseApproximation operator
      level.Request("Ainv", invapproxFact.get());

      // generate Schur complement operator
      invapproxFact->Build(level);

      RCP<Vector> Ainv = level.Get<RCP<Vector> >("Ainv", invapproxFact.get());
      TEST_EQUALITY(Ainv.is_null(), false);
      TEST_EQUALITY(Ainv->getMap()->getMinGlobalIndex(), comm->getRank() * int(n*n/comm->getSize()));
      TEST_EQUALITY(Ainv->getMap()->getMaxGlobalIndex(), comm->getRank() * int(n*n/comm->getSize()) + int(n*n/comm->getSize()-1));

      Teuchos::ArrayRCP<const Scalar> AinvData = Ainv->getData(0);
      bool bCheck = false;
      for(int i=0; i<int(n*n/comm->getSize()); i++)
        if(std::abs(AinvData[i] - Teuchos::as<Scalar>(0.166667)) < std::abs(Teuchos::as<Scalar>(1e-6)) ||
           std::abs(AinvData[i] - Teuchos::as<Scalar>(0.142857)) < std::abs(Teuchos::as<Scalar>(1e-6)))
          bCheck = true;
      TEST_EQUALITY(bCheck, true);
    }

  } //InverseLumpingConstructor

#  define MUELU_ETI_GROUP(SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InverseApproximationFactory, InverseDiagonalConstructor, SC, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InverseApproximationFactory, InverseLumpingConstructor, SC, LO, GO, Node) \

#include <MueLu_ETI_4arg.hpp>
}

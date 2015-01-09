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
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_TPETRA
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_TpetraOperator.hpp"
#endif
#ifdef HAVE_MUELU_EPETRA
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_EpetraOperator.hpp"
#endif

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  typedef MueLu::Utils<SC,LO,GO,NO> Utils;

  TEUCHOS_UNIT_TEST(PetraOperator, CreatePreconditioner)
  {
    out << "version: " << MueLu::Version() << std::endl;

    using Teuchos::RCP;

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO nx = 1000;

    // Matrix
    RCP<Matrix>     Op  = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx * comm->getSize(), lib);
    RCP<const Map > map = Op->getRowMap();

    // Normalized RHS
    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(map, 1);
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);

    // Zero initial guess
    RCP<MultiVector> X1   = MultiVectorFactory::Build(Op->getRowMap(), 1);
    X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", Op->getRowMap(), galeriList);
    RCP<MultiVector> nullspace   = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Op->getDomainMap(), 1);
    nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());
#endif

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      std::string xmlFileName = "test.xml";

      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > tpA = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Op);

      RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName);

      tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      xmlFileName = "testWithRebalance.xml";

      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpcoordinates = Utils::MV2NonConstTpetraMV(coordinates);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpnullspace   = Utils::MV2NonConstTpetraMV(nullspace);

      tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName, tpcoordinates);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

      tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName, tpcoordinates, tpnullspace);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#endif

#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Tpetra is not available");
#endif

    } else if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      std::string xmlFileName = "test.xml";

      RCP<Epetra_CrsMatrix> epA = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstEpetraCrs(Op);

      RCP<MueLu::EpetraOperator> eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName);

      eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      xmlFileName = "testWithRebalance.xml";

      RCP<Epetra_MultiVector> epcoordinates = Utils::MV2NonConstEpetraMV(coordinates);
      RCP<Epetra_MultiVector> epnullspace   = Utils::MV2NonConstEpetraMV(nullspace);

      eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName, epcoordinates);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

      eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName, epcoordinates, epnullspace);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#endif

#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Epetra is not available");
#endif

    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Unknown Xpetra lib");
    }

  } //CreatePreconditioner

  TEUCHOS_UNIT_TEST(PetraOperator, CreatePreconditioner_PDESystem)
  {
    out << "version: " << MueLu::Version() << std::endl;

    using Teuchos::RCP;

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO nx = 972;

    // Matrix
    RCP<Matrix>     Op  = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx * comm->getSize(), lib);
    RCP<const Map > map = Op->getRowMap();

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
    Teuchos::ParameterList clist;
    clist.set("nx", (nx * comm->getSize())/3);
    RCP<const Map>   cmap        = MapFactory::Build(lib, Teuchos::as<size_t>((nx * comm->getSize())/3), Teuchos::as<int>(0), comm);
    RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", cmap, clist);
    RCP<MultiVector> nullspace   = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Op->getDomainMap(), 1);
    nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());


    for (int k = 0; k < 2; k++) {
      std::string xmlFileName;
      if (k == 0) xmlFileName = "testPDE.xml";
      if (k == 1) xmlFileName = "testPDE1.xml";

      // Normalized RHS
      RCP<MultiVector> RHS1 = MultiVectorFactory::Build(Op->getRowMap(), 1);
      RHS1->setSeed(846930886);
      RHS1->randomize();
      Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
      RHS1->norm2(norms);
      RHS1->scale(1/norms[0]);

      // Zero initial guess
      RCP<MultiVector> X1   = MultiVectorFactory::Build(Op->getRowMap(), 1);
      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

      if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> >   tpA           = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Op);
        RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpcoordinates = Utils::MV2NonConstTpetraMV(coordinates);
        RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpnullspace   = Utils::MV2NonConstTpetraMV(nullspace);

        RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName, tpcoordinates);

        tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
        out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
            Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

        tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName, tpcoordinates, tpnullspace);

        X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
        tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
        out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
            Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Tpetra is not available");
#endif

      } else if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
        RCP<Epetra_CrsMatrix>   epA           = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstEpetraCrs(Op);
        RCP<Epetra_MultiVector> epcoordinates = Utils::MV2NonConstEpetraMV(coordinates);
        RCP<Epetra_MultiVector> epnullspace   = Utils::MV2NonConstEpetraMV(nullspace);

        RCP<MueLu::EpetraOperator> eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName, epcoordinates);

        eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
        out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
            Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

        eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName, epcoordinates, epnullspace);

        X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
        eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
        out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
            Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Epetra is not available");
#endif
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Unknown Xpetra lib");
      }
    }
#endif // defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
  }

  TEUCHOS_UNIT_TEST(PetraOperator, ReusePreconditioner)
  {
    out << "version: " << MueLu::Version() << std::endl;

    using Teuchos::RCP;

    Xpetra::UnderlyingLib          lib  = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    GO nx = 1000;

    // Matrix
    RCP<Matrix>     Op  = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx * comm->getSize(), lib);
    RCP<const Map > map = Op->getRowMap();

    // Normalized RHS
    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);

    // Zero initial guess
    RCP<MultiVector> X1   = MultiVectorFactory::Build(Op->getRowMap(), 1);
    X1->putScalar(Teuchos::ScalarTraits<SC>::zero());

    std::string xmlFileName = "testReuse.xml";
    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > tpA = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Op);

      RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA, xmlFileName);

      tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

      // Reuse preconditioner
      MueLu::ReuseTpetraPreconditioner(tpA, *tH);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      tH->apply(*(Utils::MV2TpetraMV(RHS1)), *(Utils::MV2NonConstTpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Tpetra is not available");
#endif

    } else if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      RCP<Epetra_CrsMatrix> epA = MueLu::Utils<SC,LO,GO,NO>::Op2NonConstEpetraCrs(Op);

      RCP<MueLu::EpetraOperator> eH = MueLu::CreateEpetraPreconditioner(epA, xmlFileName);

      eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

      // Reuse preconditioner
      MueLu::ReuseEpetraPreconditioner(epA, *eH);

      X1->putScalar(Teuchos::ScalarTraits<SC>::zero());
      eH->Apply(*(Utils::MV2EpetraMV(RHS1)), *(Utils::MV2NonConstEpetraMV(X1)));
      out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) <<
          Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Epetra is not available");
#endif

    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::InvalidArgument, "Unknown Xpetra lib");
    }
  }

}//namespace MueLuTests

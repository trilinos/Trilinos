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
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

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
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  typedef MueLu::Utils<SC,LO,GO,NO,LMO> Utils;
  typedef MueLu::TpetraOperator<SC,LO,GO,NO,LMO> TpetraOperator;

TEUCHOS_UNIT_TEST(TpetraOperator, CreatePreconditioner)
{

  out << "version: " << MueLu::Version() << std::endl;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra )
  {
    //matrix
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO nx=1000;
    RCP<Matrix> Op = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx*comm->getSize());
    RCP<const Map > map = Op->getRowMap();

    // ------------- test Tpetra Operator wrapping MueLu hierarchy ------------
    RCP< Tpetra::CrsMatrix<SC, LO, GO, NO> > tpA = MueLu::Utils<SC,LO,GO,NO,LMO>::Op2NonConstTpetraCrs(Op);
    std::string xmlFileName="test.xml";
    Teuchos::RCP<TpetraOperator> tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA,xmlFileName);

    RCP<MultiVector> RHS1 = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X1   = MultiVectorFactory::Build(Op->getRowMap(), 1);

    //normalized RHS, zero initial guess
    RHS1->setSeed(846930886);
    RHS1->randomize();
    Teuchos::Array<ST::magnitudeType> norms(1);
    RHS1->norm2(norms);
    RHS1->scale(1/norms[0]);

    X1->putScalar( (SC) 0.0);
    tH->apply(*(Utils::MV2TpetraMV(RHS1)),*(Utils::MV2NonConstTpetraMV(X1)));
    out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

    Teuchos::ParameterList galeriList;
    galeriList.set("nx", nx);
    RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",Op->getRowMap(),galeriList);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpcoordinates = Utils::MV2NonConstTpetraMV(coordinates);
    tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA,xmlFileName,tpcoordinates);

    X1->putScalar( (SC) 0.0);
    tH->apply(*(Utils::MV2TpetraMV(RHS1)),*(Utils::MV2NonConstTpetraMV(X1)));
    out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;

    RCP<Xpetra::MultiVector<SC, LO, GO, NO> > nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Op->getDomainMap(),1);
    nullspace->putScalar( Teuchos::ScalarTraits<SC>::one() );
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpnullspace = Utils::MV2NonConstTpetraMV(nullspace);
    tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(tpA,xmlFileName,tpcoordinates,tpnullspace);

    X1->putScalar( (SC) 0.0);
    tH->apply(*(Utils::MV2TpetraMV(RHS1)),*(Utils::MV2NonConstTpetraMV(X1)));
    out << "after apply, ||b-A*x||_2 = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << Utils::ResidualNorm(*Op, *X1, *RHS1) << std::endl;


  } else {

    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;

  }

} //CreatePreconditioner

}//namespace MueLuTests

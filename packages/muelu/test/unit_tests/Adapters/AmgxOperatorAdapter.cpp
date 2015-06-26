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
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_AMGXOperator.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_TpetraOperator.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

  typedef MueLu::Utils<SC,LO,GO,NO> Utils;
  typedef MueLu::AMGXOperator<double,int,int,NO> AMGXOperator;

  TEUCHOS_UNIT_TEST(AMGXOperator, Apply)
  {

    out << "version: " << MueLu::Version() << std::endl;

    if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra)
    {
      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
      if(!(comm->getSize() > 1)){
      //matrix
      int nx = 90;
      RCP<Matrix> Op = TestHelpers::TestFactory<double, int, int, NO>::Build2DPoisson(nx, -1, Xpetra::UseTpetra); 
      Teuchos::ParameterList mtxList;
      mtxList.set("nx", nx);
      mtxList.set("ny", nx);
      mtxList.set("matrixType","Laplace2D");
      //RCP<Matrix> Op = TestHelpers::TestFactory<double, int, int, NO>::BuildMatrix(mtxList, Xpetra::UseTpetra); 
      Utils::Write("A_amgx.mm", *Op);
      //RCP<const Map > map = Op->getRowMap();
      RCP<Tpetra::CrsMatrix<double, int, int,NO> > tpA = MueLu::Utils<double, int, int ,NO>::Op2NonConstTpetraCrs(Op);

      Teuchos::ParameterList params;
      params.set("use external multigrid package", "amgx");
      Teuchos::ParameterList subList = params.sublist("amgx:params", false);
      //subList.set("json file", "test.json");
      params.sublist("amgx:params").set("json file", "test.json");
       
      //subList.set("config_version", "2");
      //subList.set("monitor_residual","1");
      //subList.set("obtain_timings","1");
      //subList.set("print_grid_stats","1");
      //subList.set("exception_handling","1");
      RCP<MueLu::TpetraOperator<double, int, int, NO> > tH = MueLu::CreateTpetraPreconditioner<double, int, int, NO>(tpA, params);
      RCP<AMGXOperator> aH = Teuchos::rcp_dynamic_cast<AMGXOperator>(tH);

      TEST_EQUALITY(aH->sizeA()==nx*nx, true);

      RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRowMap(), 1);
      RCP<MultiVector> X   = MultiVectorFactory::Build(Op->getRowMap(), 1);

      //normalized RHS, zero initial guess
      //Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1); 
      //RHS->setSeed(846930886);
      //RHS->randomize();
      //RHS->norm2(norms);
      //RHS->scale(1/norms[0]);
      RHS->putScalar( (double) 1.0);
      X->putScalar( (double) 0.0);
 
      aH->apply(*(Utils::MV2TpetraMV(RHS)),*(Utils::MV2NonConstTpetraMV(X)));
      std::cout<<" status of solve: " << aH->getStatus() << " number of iterations for solve: " << aH->iters() << std::endl;
      TEST_EQUALITY(aH->getStatus()==0, true);
      }
    } else {

      out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;

    }

  } //Apply

}//namespace MueLuTests

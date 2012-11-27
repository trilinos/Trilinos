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
/*
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 */


int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor m(myTime);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // custom parameters
  LO maxLevels = 4;

  GO maxCoarseSize=1; //FIXME clp doesn't like long long int
  std::string aggOrdering = "natural";
  int minPerAgg=3;
  int maxNbrAlreadySelected=0;

  int globalNumDofs = 8898;
  int nProcs = comm->getSize();
  int nDofsPerNode = 3;

  int nLocalDofs = (int) globalNumDofs / nProcs;
  nLocalDofs = nLocalDofs - (nLocalDofs % nDofsPerNode);
  int nCumulatedDofs = 0;
  sumAll(comm,nLocalDofs, nCumulatedDofs);
  //Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM, 1, nLocalDofs, &nCumulatedDofs );

  if(comm->getRank() == nProcs-1) {
    nLocalDofs += globalNumDofs - nCumulatedDofs;
  }

  std::cout << "PROC: " << comm->getRank() << " numLocalDofs=" << nLocalDofs << std::endl;

  // read in problem
  Epetra_Map emap (globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  std::cout << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("A5932_re1000.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("b5932_re1000.txt",emap,ptrf);
  //EpetraExt::MatrixMarketFileToMultiVector( "stru2d_ns.txt", emap, ptrNS);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

  // Epetra_CrsMatrix -> Xpetra::Matrix
  RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
  RCP<Matrix> Op = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);

  Op->SetFixedBlockSize(nDofsPerNode);   // 2 velocity dofs and 1 pressure dof per node.

  // Epetra_Vector -> Xpetra::Vector
  RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVector(epNS));

  // Epetra_Map -> Xpetra::Map
  const RCP< const Map> map = Xpetra::toXpetra(emap);

  RCP<Hierarchy> H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize(maxCoarseSize);

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  //Finest->Set("Nullspace",xNS);

  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  dropFact->SetVerbLevel(MueLu::Extreme);


  //RCP<PreDropFunctionConstVal> predrop = rcp(new PreDropFunctionConstVal(0.00001));
  //dropFact->SetPreDropFunction(predrop);
  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetFactory("Graph", dropFact);
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
    *out << "aggregate ordering :                    RANDOM" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
    *out << "aggregate ordering :                    GRAPH" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  // build transfer operators
  RCP<PgPFactory> Pfact = rcp( new PgPFactory() );
  RCP<Factory> Rfact  = rcp( new GenericRFactory());
  //RCP<SaPFactory> Pfact  = rcp( new SaPFactory() );
  //RCP<Factory>   Rfact  = rcp( new TransPFactory() );
  RCP<RAPFactory> Acfact = rcp( new RAPFactory() );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 3);
  ifpackList.set("relaxation: damping factor", (SC) 0.6); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  // create coarsest smoother
  RCP<SmootherPrototype> coarsestSmooProto;
  std::string type = "";
  Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
  coarsestSmooProto = Teuchos::rcp( new DirectSolver("Superlu", coarsestSmooList) );
#else
  coarsestSmooProto = Teuchos::rcp( new DirectSolver("Klu", coarsestSmooList) );
#endif
  RCP<SmootherFactory> coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));

  FactoryManager M;
  M.SetFactory("Graph", dropFact);
  M.SetFactory("Aggregates", UCAggFact);
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarsestSmooFact);

  H->Setup(M, 0, maxLevels);

  Finest->print(*out);

  RCP<Level> coarseLevel = H->GetLevel(1);
  coarseLevel->print(*out);

  RCP<Level> coarseLevel2 = H->GetLevel(2);
  coarseLevel2->print(*out);

  RCP<MultiVector> xLsg = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    // calculate initial (absolute) residual
    Teuchos::Array<ST::magnitudeType> norms(1);
    xRhs->norm2(norms);
    *out << "||x_0|| = " << norms[0] << std::endl;

    // apply ten multigrid iterations
    H->Iterate(*xRhs,10,*xLsg);

    // calculate and print residual
    RCP<MultiVector> xTmp = MultiVectorFactory::Build(map,1);
    Op->apply(*xLsg,*xTmp,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    xRhs->update((SC)-1.0,*xTmp,(SC)1.0);
    xRhs->norm2(norms);
    *out << "||x|| = " << norms[0] << std::endl;

  }

  // TODO: don't forget to add Aztec as prerequisite in CMakeLists.txt!
  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //
  {
    RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
    X->PutScalar(0.0);
    Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

    AztecOO aztecSolver(epetraProblem);
    aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

    MueLu::EpetraOperator aztecPrec(H);
    aztecSolver.SetPrecOperator(&aztecPrec);

    int maxIts = 50;
    double tol = 1e-8;

    aztecSolver.Iterate(maxIts, tol);
  }

  return EXIT_SUCCESS;
}

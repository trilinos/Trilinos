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
 * Structure2D_epetra.cpp
 *
 *  Created on: Oct 24, 2011
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

#include <MueLu_MLParameterListInterpreter.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <MueLu_EpetraOperator.hpp>

/*!
 *  2d structural mechanics example for Epetra
 *
 *  (Nearly) Symmetric problem (except of Dirichlet boundaries) solved with AMG solver using a
 *  3 level multigrid with smoothed aggregation transfer operators.
 *
 */

void FillMLParameterList(Teuchos::ParameterList & params) {

  params.set("PDE equations",2);
  params.set("aggregation: damping factor", 1.33);
  params.set("aggregation: nodes per aggregate", 27);
  params.set("aggregation: threshold", 0.0);
  params.set("aggregation: type", "Uncoupled");
  params.set("coarse: max size", 50);
  params.set("coarse: pre or post", "post");
  params.set("coarse: sweeps", 1);
#if defined(HAVE_AMESOS_SUPERLU)
  params.set("coarse: type", "Amesos-Superlu");
#else
  params.set("coarse: type", "Amesos-KLU");
#endif
  params.set("max levels", 7);
  params.set("prec type","MGV");
  Teuchos::ParameterList & l0 = params.sublist("smoother: list (level 0)");
  Teuchos::ParameterList & l1 = params.sublist("smoother: list (level 1)");
  Teuchos::ParameterList & l2 = params.sublist("smoother: list (level 2)");
  Teuchos::ParameterList & l3 = params.sublist("smoother: list (level 3)");
  Teuchos::ParameterList & l4 = params.sublist("smoother: list (level 4)");
  Teuchos::ParameterList & l5 = params.sublist("smoother: list (level 5)");
  Teuchos::ParameterList & l6 = params.sublist("smoother: list (level 6)");


  l0.set("smoother: damping factor", 1.0);
  l0.set("smoother: sweeps", 1);
  l0.set("smoother: pre or post", "both");
  l0.set("smoother: type", "symmetric Gauss-Seidel");
  l1.set("smoother: damping factor", 1.0);
  l1.set("smoother: sweeps", 1);
  l1.set("smoother: pre or post", "both"); // "post"
  l1.set("smoother: type", "symmetric Gauss-Seidel");
  l2.set("smoother: damping factor", 1.0);
  l2.set("smoother: sweeps", 1);
  l2.set("smoother: pre or post", "both");
  l2.set("smoother: type", "symmetric Gauss-Seidel");
  l3.set("smoother: damping factor", 1.0);
  l3.set("smoother: sweeps", 1);
  l3.set("smoother: pre or post", "pre");
  l3.set("smoother: type", "symmetric Gauss-Seidel");
  l4.set("smoother: damping factor", 0.89);
  l4.set("smoother: sweeps", 1);
  l4.set("smoother: type", "Jacobi");
  l5.set("smoother: damping factor", 0.89);
  l5.set("smoother: sweeps", 12);
  l5.set("smoother: type", "Jacobi");
  l6.set("smoother: damping factor", 0.89);
  l6.set("smoother: sweeps", 14);
  l6.set("smoother: type", "Jacobi");

}

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
  Teuchos::TimeMonitor M(myTime);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif


  int globalNumDofs = 7020; //3402;
  int nProcs = comm->getSize();
  int nDofsPerNode = 2;

  int nLocalDofs = (int) globalNumDofs / nProcs;
  nLocalDofs = nLocalDofs - (nLocalDofs % nDofsPerNode);
  int nCumulatedDofs = 0;
  sumAll(comm,nLocalDofs, nCumulatedDofs);

  if(comm->getRank() == nProcs-1) {
    nLocalDofs += globalNumDofs - nCumulatedDofs;
  }

  std::cout << "PROC: " << comm->getRank() << " numLocalDofs=" << nLocalDofs << std::endl;

  // read in problem
  Epetra_Map emap (globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
  //Epetra_Map emap(3402,0,*Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  std::cout << "Reading matrix market file" << std::endl;
  /*EpetraExt::MatrixMarketFileToCrsMatrix("/home/tobias/trilinos/Trilinos_dev/ubuntu_openmpi/preCopyrightTrilinos/muelu/example/Structure/stru2d_A.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("/home/tobias/trilinos/Trilinos_dev/ubuntu_openmpi/preCopyrightTrilinos/muelu/example/Structure/stru2d_b.txt",emap,ptrf);
  EpetraExt::MatrixMarketFileToMultiVector( "/home/tobias/trilinos/Trilinos_dev/ubuntu_openmpi/preCopyrightTrilinos/muelu/example/Structure/stru2d_ns.txt", emap, ptrNS);*/
  /*EpetraExt::MatrixMarketFileToCrsMatrix("/home/wiesner/trilinos/Trilinos_dev/fc8_openmpi_dbg_q52011/preCopyrightTrilinos/muelu/example/Structure/stru2d_A.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("/home/wiesner/trilinos/Trilinos_dev/fc8_openmpi_dbg_q52011/preCopyrightTrilinos/muelu/example/Structure/stru2d_b.txt",emap,ptrf);
  EpetraExt::MatrixMarketFileToMultiVector( "/home/wiesner/trilinos/Trilinos_dev/fc8_openmpi_dbg_q52011/preCopyrightTrilinos/muelu/example/Structure/stru2d_ns.txt", emap, ptrNS);*/
  EpetraExt::MatrixMarketFileToCrsMatrix("stru2d_A.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("stru2d_b.txt",emap,ptrf);
  EpetraExt::MatrixMarketFileToMultiVector( "stru2d_ns.txt", emap, ptrNS);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);
  
  // Epetra_CrsMatrix -> Xpetra::Matrix
  RCP<CrsMatrix> exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
  RCP<Matrix> Op = crsOp;

  // Epetra_Vector -> Xpetra::Vector
  RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVector(epNS));

  // Epetra_Map -> Xpetra::Map
  const RCP< const Map> map = Xpetra::toXpetra(emap);

  Teuchos::ParameterList mlParams;
  FillMLParameterList(mlParams); // fill ML parameter list (without nullspace)

  MLParameterListInterpreter mueLuFactory(mlParams);
  RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
  H->GetLevel(0)->Set("A", Op);
  H->GetLevel(0)->Set("Nullspace", xNS);
  mueLuFactory.SetupHierarchy(*H);

  H->SetVerbLevel(MueLu::High);

  RCP<MultiVector> xLsg = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    H->Iterate(*xRhs,10,*xLsg);

    //xLsg->describe(*out,Teuchos::VERB_EXTREME);
  }

  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //
  {
    RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
    X->PutScalar(0.0);
    Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

    AztecOO aztecSolver(epetraProblem);
    aztecSolver.SetAztecOption(AZ_solver, AZ_cg);

    MueLu::EpetraOperator aztecPrec(H);
    aztecSolver.SetPrecOperator(&aztecPrec);

    int maxIts = 50;
    double tol = 1e-8;

    aztecSolver.Iterate(maxIts, tol);
  }

  /*for(int i=0; i<H->GetNumLevels(); i++) {
    RCP<Level> l = H->GetLevel(i);
    *out << std::endl << "Level " << i << std::endl;
    l->print(*out);
  }*/


  return EXIT_SUCCESS;
}

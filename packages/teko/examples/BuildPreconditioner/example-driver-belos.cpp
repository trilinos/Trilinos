// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "mpi.h"

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosThyraAdapter.hpp"  // Requires Stratimikos...

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

RCP<Teuchos::ParameterList> buildLibPL();

int main(int argc, char* argv[]) {
  typedef double ST;
  typedef Thyra::MultiVectorBase<ST> MV;
  typedef Thyra::LinearOpBase<ST> OP;

  typedef Tpetra::Vector<ST> TP_Vec;
  typedef Tpetra::CrsMatrix<ST> TP_Crs;
  typedef Tpetra::Operator<ST> TP_Op;

  typedef TP_Vec::local_ordinal_type LO;
  typedef TP_Vec::global_ordinal_type GO;
  typedef TP_Vec::node_type NT;

  // calls MPI_Init and MPI_Finalize
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder =
      Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);

  // Build the Tpetra matrices and vectors
  /////////////////////////////////////////////////////////

  // read in the CRS matrix
  RCP<TP_Crs> crsMat = Tpetra::MatrixMarket::Reader<TP_Crs>::readSparseFile(
      "../data/nsjac.mm", Tpetra::getDefaultComm());

  RCP<TP_Crs> zeroCrsMat = rcp(new TP_Crs(*crsMat, Teuchos::Copy));
  zeroCrsMat->setAllToScalar(0.0);

  RCP<TP_Op> Mat     = crsMat;
  RCP<TP_Op> zeroMat = zeroCrsMat;

  // Allocate some right handside vectors
  RCP<TP_Vec> x0_tp = rcp(new TP_Vec(Mat->getDomainMap()));
  RCP<TP_Vec> x1_tp = rcp(new TP_Vec(Mat->getDomainMap()));
  RCP<TP_Vec> b0_tp = rcp(new TP_Vec(Mat->getRangeMap()));
  RCP<TP_Vec> b1_tp = rcp(new TP_Vec(Mat->getRangeMap()));
  b0_tp->randomize();
  b1_tp->randomize();

  RCP<const Thyra::TpetraVectorSpace<ST, LO, GO, NT> > domain =
      Thyra::tpetraVectorSpace<ST>(Mat->getDomainMap());
  RCP<const Thyra::TpetraVectorSpace<ST, LO, GO, NT> > range =
      Thyra::tpetraVectorSpace<ST>(Mat->getRangeMap());

  // Build Teko compatible matrices and vectors
  /////////////////////////////////////////////////////////

  // convert them to teko compatible sub vectors
  Teko::MultiVector x0_th = Thyra::tpetraVector(domain, x0_tp);
  Teko::MultiVector x1_th = Thyra::tpetraVector(domain, x1_tp);
  Teko::MultiVector b0_th = Thyra::tpetraVector(range, b0_tp);
  Teko::MultiVector b1_th = Thyra::tpetraVector(range, b1_tp);
  std::vector<Teko::MultiVector> x_vec;
  x_vec.push_back(x0_th);
  x_vec.push_back(x1_th);
  std::vector<Teko::MultiVector> b_vec;
  b_vec.push_back(b0_th);
  b_vec.push_back(b1_th);

  Teko::MultiVector x =
      Teko::buildBlockedMultiVector(x_vec);  // these will be used in the Teko solve
  Teko::MultiVector b = Teko::buildBlockedMultiVector(b_vec);

  // Build the Teko compatible linear system
  Teko::LinearOp thMat  = Thyra::tpetraLinearOp<double>(range, domain, Mat);
  Teko::LinearOp thZero = Thyra::tpetraLinearOp<double>(range, domain, zeroMat);
  Teko::LinearOp A =
      Thyra::block2x2(thMat, thZero, thZero, thMat);  // build an upper triangular 2x2

  // Build the preconditioner
  /////////////////////////////////////////////////////////

  // build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib =
      Teko::InverseLibrary::buildFromParameterList(*buildLibPL(), linearSolverBuilder);

  // build the inverse factory needed by the example preconditioner
  RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("Gauss-Seidel");

  // build the preconditioner from the jacobian
  Teko::LinearOp prec = Teko::buildInverse(*inverse, A);

  // Setup the Belos solver
  /////////////////////////////////////////////////////////

  Teuchos::ParameterList belosList;
  belosList.set("Num Blocks", 200);              // Maximum number of blocks in Krylov factorization
  belosList.set("Block Size", 1);                // Blocksize to be used by iterative solver
  belosList.set("Maximum Iterations", 200);      // Maximum number of iterations allowed
  belosList.set("Maximum Restarts", 1);          // Maximum number of restarts allowed
  belosList.set("Convergence Tolerance", 1e-5);  // Relative convergence tolerance requested
  belosList.set(
      "Verbosity",
      33);  // Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
  belosList.set("Output Frequency", 1);
  belosList.set("Output Style", 1);

  RCP<Belos::LinearProblem<double, MV, OP> > problem =
      rcp(new Belos::LinearProblem<double, MV, OP>(A, x, b));
  problem->setLeftPrec(prec);
  problem->setProblem();  // should check the return type!!!

  RCP<Belos::SolverManager<double, MV, OP> > solver =
      rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcpFromRef(belosList)));

  //
  // Perform solve
  //
  Belos::ReturnType ret = solver->solve();

  if (ret != Belos::Converged) {
    std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    return -1;
  }

  return 0;
}

RCP<Teuchos::ParameterList> buildLibPL() {
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  {
    Teuchos::ParameterList& sub_jac = pl->sublist("Jacobi");
    sub_jac.set("Type", "Block Jacobi");
    sub_jac.set("Inverse Type", "Ifpack2");

    Teuchos::ParameterList& sub_gs = pl->sublist("Gauss-Seidel");
    sub_gs.set("Type", "Block Gauss-Seidel");
    sub_gs.set("Use Upper Triangle", true);
    sub_gs.set("Inverse Type", "Ifpack2");
  }
  return pl;
}

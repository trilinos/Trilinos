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
#include <iostream>

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Epetra includes
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

// // EpetraExt includes
// #include <EpetraExt_CrsMatrixIn.h>
// #include <EpetraExt_VectorIn.h>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// MueLu includes
#include <MueLu.hpp> // TODO Usefull?
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

int main(int argc,char * argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization
  //

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // build global communicator TODO: convert from Teuchos::Comm
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false);

  Galeri::Xpetra::Parameters<int> matrixParameters(clp, 256); // manage parameters of the test case
  // Xpetra::Parameters              xpetraParameters(clp);   // manage parameters of xpetra
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra; // Epetra only for the moment

  std::string xmlFileName = "stratimikos_ParameterList.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'stratimikos_ParameterList.xml'.");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  // Read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);

  //
  // Construct the problem
  //

  //   // Read in the matrix, store pointer as an RCP
  //   Epetra_CrsMatrix * ptrA = 0;
  //   EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac_test.mm",Comm,ptrA);
  //   RCP<Epetra_CrsMatrix> A = rcp(ptrA);
  //
  //   // read in the RHS vector
  //   Epetra_Vector * ptrb = 0;
  //   EpetraExt::MatrixMarketFileToVector("../data/nsrhs_test.mm",A->MatrixRangeMap(),ptrb);
  //   RCP<const Epetra_Vector> b = rcp(ptrb);

  RCP<const Map> map = MapFactory::createUniformContigMap(lib, matrixParameters.GetNumGlobalElements(), comm);
  RCP<Galeri::Xpetra::Problem<Map,Xpetra::EpetraCrsMatrix,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<double, int, int, Map,  Xpetra::EpetraCrsMatrix, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
  RCP<const Epetra_CrsMatrix> A = Pr->BuildMatrix()->getEpetra_CrsMatrix();

  //
  // Allocate vectors
  //

  RCP<Epetra_Vector> X = rcp(new Epetra_Vector(A->DomainMap()));
  X->PutScalar(0.0);

  RCP<Epetra_Vector> B = rcp(new Epetra_Vector(A->DomainMap()));
  B->SetSeed(846930886); B->Random();

  //
  // Build Thyra linear algebra objects
  //

  RCP<const Thyra::LinearOpBase<double> > thyraA = Thyra::epetraLinearOp(A);
  RCP<const Thyra::VectorBase<double> >   thyraB = Thyra::create_Vector(B, thyraA->range());
  RCP<Thyra::VectorBase<double> >         thyraX = Thyra::create_Vector(X, thyraA->domain());

  //
  // Build Stratimikos solver
  //

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
  Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder);     // Register MueLu as a Stratimikos preconditioner strategy.
  linearSolverBuilder.setParameterList(paramList);              // Setup solver parameters using a Stratimikos parameter list.

  // Build a new "solver factory" according to the previously specified parameter list.
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

  // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > thyraInverseA = Thyra::linearOpWithSolve(*solverFactory, thyraA);

  //
  // Solve Ax = b.
  //

  Thyra::SolveStatus<double> status = Thyra::solve<double>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
  std::cout << status << std::endl;

  return (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED) ? EXIT_SUCCESS : EXIT_FAILURE;
}

// Thyra::assign(thyraX.ptr(), 0.0);

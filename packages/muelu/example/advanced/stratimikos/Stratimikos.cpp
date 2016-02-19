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
#include <iostream>

/*
   Call MueLu via the Stratimikos interface.

Usage:
./MueLu_Stratimikos.exe : use xml configuration file stratimikos_ParameterList.xml

Note:
The source code is not MueLu specific and can be used with any Stratimikos strategy.
*/

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra includes
#include <Epetra_Vector.h>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// MueLu includes
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <MueLu_UseDefaultTypes.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

int main(int argc,char * argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization
  //

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    //
    // Parameters
    //

    Teuchos::CommandLineProcessor clp(false);

    Galeri::Xpetra::Parameters<int> matrixParameters(clp, 256); // manage parameters of the test case
    Xpetra::Parameters              xpetraParameters(clp);   // manage parameters of xpetra

    std::string xmlFileName = "stratimikos_ParameterList.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'stratimikos_ParameterList.xml'.");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    // Read in parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);

    //
    // Construct the problem
    //
    RCP<const Map> map = MapFactory::createUniformContigMap(lib, matrixParameters.GetNumGlobalElements(), comm);

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<CrsMatrixWrap> A = Pr->BuildMatrix();

    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);

    {
      // we set seed for reproducibility
      Utilities::SetRandomSeed(*comm);
      X->randomize();
      A->apply(*X, *B, Teuchos::NO_TRANS, Teuchos::ScalarTraits<Scalar>::one(), Teuchos::ScalarTraits<Scalar>::zero());

      Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(Teuchos::ScalarTraits<Scalar>::one()/norms[0]);
      X->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    }

    //
    // Build Thyra linear algebra objects
    //

    RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(A->getCrsMatrix());

    RCP<const Thyra::VectorBase<Scalar> > thyraB = Teuchos::null;
    RCP<      Thyra::VectorBase<Scalar> > thyraX = Teuchos::null;
#ifdef HAVE_XPETRA_TPETRA
    if(lib==Xpetra::UseTpetra) {
      thyraB = Thyra::createVector(Xpetra::toTpetra(*B), thyraA->range());
      thyraX = Thyra::createVector(Xpetra::toTpetra(*X), thyraA->domain());
    }
#endif
#ifdef HAVE_XPETRA_EPETRA
    if(lib==Xpetra::UseEpetra) {
      RCP<Epetra_Vector> epVecB = Teuchos::rcpFromRef<Epetra_Vector>(Xpetra::toEpetra(*B));
      RCP<Epetra_Vector> epVecX = Teuchos::rcpFromRef<Epetra_Vector>(Xpetra::toEpetra(*X));
      thyraB = Thyra::create_Vector(epVecB, thyraA->range());
      thyraX = Thyra::create_Vector(epVecX, thyraA->domain());
    }
#endif


    //
    // Build Stratimikos solver
    //

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::enableMueLu<LocalOrdinal,GlobalOrdinal,Node>(linearSolverBuilder);                // Register MueLu as a Stratimikos preconditioner strategy.
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

    success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

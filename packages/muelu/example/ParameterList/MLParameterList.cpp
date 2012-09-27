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

#include <MueLu_ConfigDefs.hpp>

#include <Teuchos_XMLParameterListHelpers.hpp> // getParametersFromXmlFile()
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
#include <Epetra_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

#ifdef HAVE_MUELU_AZTECOO
#include <AztecOO.h>
#endif

#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#endif

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

// Default problem is Laplace1D with nx = 8748. Use --help to list available options.

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  //TODO: FIXME: option by default does not work for MueLu/Tpetra

  int nIts = 9;

  Teuchos::CommandLineProcessor clp(false); // Note: 

  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 256); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);      // manage parameters of xpetra

  std::string xmlFileName; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default an hard-coded parameter list.");
  int muelu = true;        clp.setOption("muelu", &muelu,       "use muelu"); //TODO: bool instead of int
  int ml    = false;       
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
  clp.setOption("ml",    &ml,          "use ml");
#endif

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  // TODO: check -ml and --linAlgebra

  if (comm->getRank() == 0) { std::cout << xpetraParameters << matrixParameters; }
  if (ml && xpetraParameters.GetLib() == Xpetra::UseTpetra) {
    ml = 0;
    std::cout << "ML preconditionner can only be built if --linAlgebra=0 (Epetra). Option --ml ignored" << std::endl;
  }

  //
  // Construct the problem
  //

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Matrix>  A   = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  //
  // Preconditionner configuration
  //
  
  // ML parameter list
  RCP<Teuchos::ParameterList> params;
  if (xmlFileName != "") {

    std::cout << "Reading " << xmlFileName << " ..." << std::endl;
    params = Teuchos::getParametersFromXmlFile(xmlFileName);

  } else {

    std::cout << "Using hard-coded parameter list:" << std::endl;
    params = rcp(new Teuchos::ParameterList());
    
    params->set("max levels", 2); //FIXME: does not work?
    
    //    params->set("smoother: damping factor", 0.9);
    // params->set("smoother: sweeps", 1);
    // params->set("smoother: pre or post", "both");
    params->set("smoother: type", "symmetric Gauss-Seidel");
    params->set("coarse: type", "Amesos-KLU");
    
    /*      Teuchos::ParameterList & l0 = params->sublist("smoother: list (level 0)");
            l0.set("smoother: damping factor", 0.9);
            l0.set("smoother: sweeps", 1);
            l0.set("smoother: pre or post", "both");
            l0.set("smoother: type", "symmetric Gauss-Seidel");
    */
    
    params->set("coarse: type","Amesos-KLU");
  }

  std::cout << *params << std::endl;
    
  if (muelu) {

    //
    // Construct a multigrid preconditioner
    //
    
    // Multigrid Hierarchy
    MLParameterListInterpreter mueLuFactory(*params);
    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    
    mueLuFactory.SetupHierarchy(*H);

    //
    // Solve Ax = b
    //
    
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);
    
    X->putScalar((Scalar) 0.0);
    B->setSeed(846930886); B->randomize();
    
    // AMG as a standalone solver
    H->IsPreconditioner(false);
    H->Iterate(*B, nIts, *X);
    
    // Print relative residual norm
    ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *X, *B)[0];
    if (comm->getRank() == 0)
      std::cout << "||Residual|| = " << residualNorms << std::endl;

#ifdef HAVE_MUELU_AZTECOO

    // AMG as a preconditioner

    //TODO: name mueluPrec and mlPrec not 

    H->IsPreconditioner(true);
    MueLu::EpetraOperator mueluPrec(H); // Wrap MueLu preconditioner into an Epetra Operator

    //
    // Solve Ax = b
    //
    RCP<Epetra_CrsMatrix> eA; //duplicate code
    { // TODO: simplify this
      RCP<CrsMatrixWrap>       xCrsOp  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
      RCP<CrsMatrix>         xCrsMtx = xCrsOp->getCrsMatrix();
      RCP<EpetraCrsMatrix>   eCrsMtx = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xCrsMtx, true);
      eA = eCrsMtx->getEpetra_CrsMatrixNonConst();
    }

    RCP<Epetra_Vector> eX = rcp(new Epetra_Vector(eA->RowMap()));
    RCP<Epetra_Vector> eB = rcp(new Epetra_Vector(eA->RowMap()));
    
    eX->PutScalar((Scalar) 0.0);
    eB->SetSeed(846930886); eB->Random();

    Epetra_LinearProblem eProblem(eA.get(), eX.get(), eB.get());

    // AMG as a standalone solver
    AztecOO solver(eProblem);
    solver.SetPrecOperator(&mueluPrec);
    solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
    solver.SetAztecOption(AZ_output, 1);
    
    solver.Iterate(nIts, 1e-10);

    { //TODO: simplify this
      RCP<Vector> mueluX = rcp(new Xpetra::EpetraVector(eX)); 
      RCP<Vector> mueluB = rcp(new Xpetra::EpetraVector(eB));
      // Print relative residual norm
      ST::magnitudeType residualNorms2 = Utils::ResidualNorm(*A, *mueluX, *mueluB)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms2 << std::endl;
    }

    // TODO: AMG as a preconditioner (AZ_cg)
#endif // HAVE_MUELU_AZTECOO

  }

#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
  if (ml) {

    std::cout << std::endl << std::endl << std::endl << std::endl << "**** ML ml ML ml ML" << std::endl << std::endl << std::endl << std::endl;

    //
    // Construct a multigrid preconditioner
    //

    // Multigrid Hierarchy
    RCP<CrsMatrixWrap>      crsOp         = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
    RCP<CrsMatrix>        crsMtx        = crsOp->getCrsMatrix();
    RCP<EpetraCrsMatrix>  epetraCrsMtx  = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(crsMtx, true);
    RCP<const Epetra_CrsMatrix> epetra_CrsMtx = epetraCrsMtx->getEpetra_CrsMatrix();

    RCP<Epetra_CrsMatrix> eA;
    { // TODO: simplify this
      RCP<CrsMatrixWrap>       xCrsOp  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
      RCP<CrsMatrix>         xCrsMtx = xCrsOp->getCrsMatrix();
      RCP<EpetraCrsMatrix>   eCrsMtx = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xCrsMtx, true);
      eA = eCrsMtx->getEpetra_CrsMatrixNonConst();
    }

    RCP<ML_Epetra::MultiLevelPreconditioner> mlPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*eA, *params));
    
#ifdef HAVE_MUELU_AZTECOO

    //
    // Solve Ax = b
    //
    
    RCP<Epetra_Vector> eX = rcp(new Epetra_Vector(eA->RowMap()));
    RCP<Epetra_Vector> eB = rcp(new Epetra_Vector(eA->RowMap()));
    
    eX->PutScalar((Scalar) 0.0);
    eB->SetSeed(846930886); eB->Random();

    Epetra_LinearProblem eProblem(eA.get(), eX.get(), eB.get());

    // AMG as a standalone solver
    AztecOO solver(eProblem);
    solver.SetPrecOperator(mlPrec.get());
    solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
    solver.SetAztecOption(AZ_output, 1);
    
    solver.Iterate(nIts, 1e-10);

    { //TODO: simplify this
      RCP<Vector> mueluX = rcp(new Xpetra::EpetraVector(eX)); 
      RCP<Vector> mueluB = rcp(new Xpetra::EpetraVector(eB));
      // Print relative residual norm
      ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *mueluX, *mueluB)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms << std::endl;
    }

    // TODO: AMG as a preconditioner (AZ_cg)
#else
    std::cout << "Enable AztecOO to see solution" << std::endl;
#endif // HAVE_MUELU_AZTECOO
  } // if (ml)

#endif // HAVE_MUELU_ML && HAVE_MUELU_EPETRA

  return EXIT_SUCCESS;
}


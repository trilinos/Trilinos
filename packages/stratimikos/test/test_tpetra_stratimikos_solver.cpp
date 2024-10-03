// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_tpetra_stratimikos_solver.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Teuchos_ParameterList.hpp"

bool Thyra::test_tpetra_stratimikos_solver(
  Teuchos::ParameterList *paramList_inout
  ,const bool dumpAll
  ,Teuchos::FancyOStream *out
  )
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::getParameter;
  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;

  bool success = true, result = false;

  try {

    TEUCHOS_TEST_FOR_EXCEPT(!paramList_inout);

    RCP<ParameterList>
      paramList = rcp(paramList_inout,false);

    if(out) {
      *out << "\nEchoing input parameters ...\n";
      paramList->print(*out,1,true,false);
    }

    // Create list of valid parameter sublists
    Teuchos::ParameterList validParamList("test_tpetra_stratimikos_solver");
    validParamList.set("Matrix File","fileName");
    validParamList.set("Solve Adjoint",false);
    validParamList.sublist("Linear Solver Builder").disableRecursiveValidation();
    validParamList.sublist("LinearOpWithSolveTester").disableRecursiveValidation();

    if(out) *out << "\nValidating top-level input parameters ...\n";
    paramList->validateParametersAndSetDefaults(validParamList);

    const std::string
      &matrixFile = getParameter<std::string>(*paramList,"Matrix File");
    const bool
      solveAdjoint = getParameter<bool>(*paramList,"Solve Adjoint");
    RCP<ParameterList>
      solverBuilderSL  = sublist(paramList,"Linear Solver Builder",true),
      lowsTesterSL     = sublist(paramList,"LinearOpWithSolveTester",true);

    if(out) *out << "\nReading in an tpetra matrix A from the file \'"<<matrixFile<<"\' ...\n";

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    RCP<const Tpetra::CrsMatrix<Scalar, LO, GO, NO>> tpetra_A = Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, NO>>::readSparseFile( matrixFile, comm );

    RCP<const LinearOpBase<double> >
      A = Thyra::createConstLinearOp(Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar,LO,GO,NO>>(tpetra_A));

    if(out) *out << "\nCreating a Stratimikos::DefaultLinearSolverBuilder object ...\n";

    RCP<Thyra::LinearSolverBuilderBase<double> >
      linearSolverBuilder = rcp(new Stratimikos::DefaultLinearSolverBuilder);

    if(out) {
      *out << "\nValid parameters for DefaultLinearSolverBuilder ...\n";
      linearSolverBuilder->getValidParameters()->print(*out,1,true,false);
    }

    linearSolverBuilder->setParameterList(solverBuilderSL);

    if(out) *out << "\nCreating the LinearOpWithSolveFactoryBase object lowsFactory ...\n";
    RCP<LinearOpWithSolveFactoryBase<double> >
      lowsFactory = createLinearSolveStrategy(*linearSolverBuilder);
    if(out) *out << "\nlowsFactory described as:\n" << describe(*lowsFactory,Teuchos::VERB_MEDIUM) << std::endl;

    if(out) *out << "\nRunning example use cases for not externally preconditioned ...\n";

    TEUCHOS_ASSERT(out != NULL);
    nonExternallyPreconditionedLinearSolveUseCases(
      *A, *lowsFactory, solveAdjoint, *out
      );

    Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
    linearOpWithSolveTester.setParameterList(lowsTesterSL);
    linearOpWithSolveTester.turn_off_all_tests();
    linearOpWithSolveTester.check_forward_default(true);
    linearOpWithSolveTester.check_forward_residual(true);
    if (solveAdjoint) {
      linearOpWithSolveTester.check_adjoint_default(true);
      linearOpWithSolveTester.check_adjoint_residual(true);
    }
    // ToDo: Use parameter lists for the above

    if(out) *out << "\nChecking the LOWSB interface ...\n";
    RCP<Thyra::LinearOpWithSolveBase<Scalar> >
      lowsA = Thyra::linearOpWithSolve<Scalar>(*lowsFactory, A);
    result = linearOpWithSolveTester.check(*lowsA, out);
    if (!result) success = false;

    if(out) {
      *out << "\nPrinting the parameter list (showing what was used) ...\n";
      paramList->print(*out,1,true,true);
    }

  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }

  return success;

}

// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "test_single_stratimikos_solver.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

bool Thyra::test_single_stratimikos_solver(
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

  bool success = true, result = false;

  try {

    TEST_FOR_EXCEPT(!paramList_inout);

    RCP<ParameterList>
      paramList = rcp(paramList_inout,false);

    if(out) {
      *out << "\nEchoing input parameters ...\n";
      paramList->print(*out,1,true,false);
    }

    // Create list of valid parameter sublists
    Teuchos::ParameterList validParamList("test_single_stratimikos_solver");
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

    if(out) *out << "\nReading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    RCP<Epetra_CrsMatrix> epetra_A;
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

    RCP<const LinearOpBase<double> >
      A = Thyra::epetraLinearOp(epetra_A);

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

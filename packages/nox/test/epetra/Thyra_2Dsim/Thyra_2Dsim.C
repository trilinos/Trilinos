//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
          
// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "EpetraModelEval2DSim.H"

using namespace std;

int main(int argc, char *argv[])
{
 
  int status = 0;
  
  // Parse the command line
  using Teuchos::CommandLineProcessor;
  CommandLineProcessor  clp;
  clp.throwExceptions(false);
  clp.addOutputSetupOptions(true);
  bool verbose = false;
  clp.setOption( "v", "disable-verbosity", &verbose, "Enable verbosity" );
  
  CommandLineProcessor::EParseCommandLineReturn
    parse_return = clp.parse(argc,argv,&std::cerr);
  
  if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
    return parse_return;

  if (verbose) 
    cout << "Verbosity Activated" << endl;
  else
    cout << "Verbosity Disabled" << endl;


  // Start up MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Check we have only one processor since this problem doesn't work
  // for more than one proc
  if (Comm.NumProc() > 1) {
    std::cerr << "Error!  Problem can only be run with at most 1 processor!"
	      << std::endl;
    return -1;
  }

  // Create the EpetraExt model evaluator object
  double d = 10.0;
  double p0 = 2.0;
  double p1 = 0.0;
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<EpetraExt::ModelEvaluator> epetraModel = 
    rcp(new EpetraModelEval2DSim(Teuchos::rcp(&Comm,false),
				 d,p0,p1,x00,x01));

  // Create the linear solver type with Stratimikos
  //Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
  //lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());

  ::Stratimikos::DefaultLinearSolverBuilder builder;
  
  Teuchos::RCP<Teuchos::ParameterList> p = 
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->set("Preconditioner Type", "Ifpack");
  //p->set("Enable Delayed Solver Construction", true);
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> > 
    lowsFactory = builder.createLinearSolveStrategy("");

  // Create the Thyra model evalutor (form the epetraext model and
  // linear solver)
  Teuchos::RCP< ::Thyra::EpetraModelEvaluator>
    epetraThyraModel = rcp(new ::Thyra::EpetraModelEvaluator());
  epetraThyraModel->initialize(epetraModel,lowsFactory);
  Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel = 
    epetraThyraModel;

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = thyraModel->getNominalValues().get_x()->clone_v();

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group = 
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel));

//   nox_group->computeF();
//   cout << "ComputedF!" << endl;
//   const NOX::Thyra::Vector& t_vec = 
//     dynamic_cast<const NOX::Thyra::Vector&>(nox_group->getF());
  
//   t_vec.print(std::cout);
//   exit(0);


  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  if (solvStatus == NOX::StatusTest::Converged)
    std::cout << "Test passed!" << std::endl;

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}

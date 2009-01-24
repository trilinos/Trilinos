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
#include "AztecOO.h"

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

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
 
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

  ::Stratimikos::DefaultLinearSolverBuilder builder;
  
  Teuchos::RCP<Teuchos::ParameterList> p = 
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->set("Preconditioner Type", "Ifpack");
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

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");

  // Create a Thyra nonlinear solver
  Teuchos::RCP< ::Thyra::NonlinearSolverBase<double> > solver = 
    Teuchos::rcp(new ::Thyra::NOXNonlinearSolver);
  
  solver->setParameterList(nl_params);
  solver->setModel(thyraModel);

  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = thyraModel->getNominalValues().get_x()->clone_v();

  ::Thyra::SolveCriteria<double> solve_criteria;
  ::Thyra::SolveStatus<double> solve_status;

  solve_status = solver->solve(initial_guess.get(), &solve_criteria);

  if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED)
    std::cout << "Test passed!" << std::endl;

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}

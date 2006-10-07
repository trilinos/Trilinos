// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <iostream>

#include "NOX.H"
#include "NOX_Epetra.H"

#include "FEApp_ModelEvaluator.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[]) {
  unsigned int nelem = 100;
  double h = 1.0/nelem;
  double alpha = 0.25;
  double beta = 0.0;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
    Teuchos::RefCountPtr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    int MyPID = Comm->MyPID();
    
    // Create mesh
    vector<double> x(nelem+1);
    for (unsigned int i=0; i<=nelem; i++)
      x[i] = h*i;

    // Set up application parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> appParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
    problemParams.set("Name", "Brusselator");
    problemParams.set("alpha", alpha);
    problemParams.set("beta", beta);
    problemParams.set("D1", D1);
    problemParams.set("D2", D2);

    // Create application
    Teuchos::RefCountPtr<FEApp::Application> app = 
      Teuchos::rcp(new FEApp::Application(x, Comm, appParams));

    // Create initial guess
    Teuchos::RefCountPtr<const Epetra_Vector> u = app->getInitialSolution();
    NOX::Epetra::Vector nox_u(*u);

    // Create model evaluator
    Teuchos::RefCountPtr<FEApp::ModelEvaluator> model = 
      Teuchos::rcp(new FEApp::ModelEvaluator(app));

    //app->getJacobianGraph()->Print(cout);

    // Create NOX interface
    Teuchos::RefCountPtr<NOX::Epetra::ModelEvaluatorInterface> interface =
      Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));

    // Set up NOX parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> noxParams =
      Teuchos::rcp(&(appParams->sublist("NOX")),false);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
		    NOX::Utils::OuterIteration + 
		    NOX::Utils::OuterIterationStatusTest + 
		    NOX::Utils::InnerIteration +
		    NOX::Utils::Parameters + 
		    NOX::Utils::Details + 
		    NOX::Utils::LinearSolverDetails +
		    NOX::Utils::Warning + 
		    NOX::Utils::Error);

    // Create printing utilities
    NOX::Utils utils(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 800);  
    lsParams.set("Tolerance", 1e-4); 
    lsParams.set("Output Frequency", 50);
    lsParams.set("Preconditioner", "Ifpack");

    // Create the Jacobian matrix
    Teuchos::RefCountPtr<Epetra_Operator> A = model->create_W(); 

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, 
							lsParams,
							iReq, iJac, A, nox_u));

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, nox_u, linsys)); 

    // Create the Solver convergence test
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> wrms = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					   NOX::StatusTest::NormF::Unscaled));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the solver
    NOX::Solver::Manager solver(grp, combo, noxParams);

    // Solve
    NOX::StatusTest::StatusType status = solver.solve();
    if (status != NOX::StatusTest::Converged)
      if (MyPID==0) 
	utils.out() << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // Output the parameter list
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << endl << "Final Parameters" << endl
		  << "****************" << endl;
      appParams->print(utils.out());
      utils.out() << endl;
    }

    //finalSolution.Print(utils.out());

#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

}

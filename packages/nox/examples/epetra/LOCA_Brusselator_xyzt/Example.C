//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

// 1D Finite Element Brusselator Test Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T    
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2   
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C    
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2   
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 *
 * This problem is examined with a variety of time integration schemes in:
 * "Studies on the Convergence of Various Time-Integration Schemes for the
 * Radiation-Diffusion Problem," Curtis C. Ober & John N. Shadid, in prep.
 *
 * In this example, only a 1st-order fully implicit (backward Euler)
 * time integration scheme is considered currently.
 *
 * Values for time step size and finite spatial extent are specified in
 * the constructor initialization list in Brusselator.C using
 * variables dt  and xmin,xmax, respectively.
 * The number of time steps to be taken is specified by variable
 * maxTimeSteps below.
 */

// NOX Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_Common.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Teuchos_GlobalMPISession.hpp"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring 
#include <vector> 
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h" 
#endif

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "Brusselator.H"              

#ifdef HAVE_NOX_EPETRAEXT
#ifdef HAVE_MPI
// Comment out following line for usual implicit time stepping on all procs
#define DO_XYZT 1
#define DO_XYZT_PREC 1
#include "EpetraExt_MultiMpiComm.h"
#endif
#endif

#ifdef DO_XYZT
#include "LOCA_Epetra_Interface_xyzt.H"              
#endif

using namespace std;

int main(int argc, char *argv[])
{
  // Initialize MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  // Get the number of elements from the command line
  int NumGlobalNodes = 100 + 1;

#ifdef DO_XYZT
  // MPI MANIPULATION FOR XYZT PROBLEMS
  int spatialProcs = 1; // default
  if (argc>2) { spatialProcs = atoi(argv[2]);}
  int numTimeSteps= 1; // default
  if (argc>3) { numTimeSteps = atoi(argv[3]);}

  Teuchos::RCP<EpetraExt::MultiMpiComm> globalComm = 
    Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD, spatialProcs, numTimeSteps));
  Epetra_Comm& Comm = globalComm->SubDomainComm();

#else
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm();
#endif

#endif

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalNodes < NumProc) {
    std::cout << "numGlobalNodes = " << NumGlobalNodes 
	 << " cannot be < number of processors = " << NumProc << std::endl;
    exit(1);
  }

  // Create the Brusselator problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (F) and Jacobian evaluation routines.
  Brusselator::OverlapType OType = Brusselator::ELEMENTS;
  Brusselator Problem(NumGlobalNodes, Comm, OType);
  double dt = 0.5;

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();
  Epetra_Vector& initCond = soln;

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list

  Teuchos::RCP<Teuchos::ParameterList> paramList =
       Teuchos::rcp(new Teuchos::ParameterList);

  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Method", "Natural");
  //locaStepperList.set("Continuation Method", "Arc Length");
  //locaStepperList.set("Continuation Method", "Householder Arc Length");
  locaStepperList.set("Continuation Parameter", "alpha");
  locaStepperList.set("Initial Value", 0.6);
  locaStepperList.set("Max Value", 100.0);
  locaStepperList.set("Min Value", 0.05);
#ifdef DO_XYZT
  locaStepperList.set("Max Steps", 7);
#else
  locaStepperList.set("Max Steps", 0);// must be 0 so just a nonlinear solver
#endif
  locaStepperList.set("Max Nonlinear Iterations", 15);
  locaStepperList.set("Enable Arc Length Scaling", true);
  locaStepperList.set("Goal Arc Length Parameter Contribution", 0.5);
  locaStepperList.set("Max Arc Length Parameter Contribution", 0.7);
  locaStepperList.set("Initial Scale Factor", 1.0);
  locaStepperList.set("Min Scale Factor", 1.0e-8);
  locaStepperList.set("Enable Tangent Factor Step Size Scaling",false);
  locaStepperList.set("Min Tangent Factor", 0.8);
  locaStepperList.set("Tangent Factor Exponent",1.5);

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
  stepSizeList.set("Method", "Constant");
 // stepSizeList.set("Method", "Adaptive");
  stepSizeList.set("Initial Step Size", -0.1);
  stepSizeList.set("Min Step Size", 1.0e-3);
  stepSizeList.set("Max Step Size", 2000.0);
  stepSizeList.set("Aggressiveness", 0.1);
  stepSizeList.set("Failed Step Reduction Factor", 0.5);
  stepSizeList.set("Successful Step Increase Factor", 1.00); // for constant

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.set("Method", "Constant");
  //predictorList.set("Method", "Tangent");
  predictorList.set("Method", "Secant");

  // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");

  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.set("Compute Eigenvalues",false);
#ifdef HAVE_LOCA_ANASAZI
  Teuchos::ParameterList& aList = locaStepperList.sublist("Eigensolver");
  aList.set("Method", "Anasazi");
  aList.set("Block Size", 1);
  aList.set("Num Blocks", 10);
  aList.set("Num Eigenvalues", 3);
  aList.set("Convergence Tolerance", 2.0e-7);
  aList.set("Convergence Check", 1);
  aList.set("Maximum Restarts",2);
  aList.set("Step Size",1);
  aList.set("Verbosity",
	    Anasazi::Errors + 
	    Anasazi::Warnings +
	    Anasazi::FinalSummary);
#endif

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			   NOX::Utils::OuterIteration + 
			   NOX::Utils::OuterIterationStatusTest + 
			   NOX::Utils::InnerIteration +
			   NOX::Utils::Parameters + 
			   NOX::Utils::Details + 
			   NOX::Utils::Warning + 
			   NOX::Utils::StepperIteration +
			   NOX::Utils::StepperDetails);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
  //searchParams.set("Method", "Interval Halving");
  //searchParams.set("Method", "Polynomial");
  //searchParams.set("Method", "NonlinearCG");
  //searchParams.set("Method", "Quadratic");
  //searchParams.set("Method", "More'-Thuente");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");
    //newtonParams.set("Forcing Term Method", "Type 1");
    //newtonParams.set("Forcing Term Method", "Type 2");
    //newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.set("Forcing Term Maximum Tolerance", 0.1);
  //dirParams.set("Method", "Steepest Descent");
  //Teuchos::ParameterList& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.set("Scaling Type", "None");
    //sdParams.set("Scaling Type", "2-Norm");
    //sdParams.set("Scaling Type", "Quadratic Model Min");
  //dirParams.set("Method", "NonlinearCG");
  //Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.set("Restart Frequency", 2000);
    //nlcgParams.set("Precondition", "On");
    //nlcgParams.set("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 800);  
  lsParams.set("Tolerance", 1e-6);
  lsParams.set("Output Frequency", 50);    
#ifdef DO_XYZT_PREC
  lsParams.set("Preconditioner", "User Defined"); 
#else
  lsParams.set("Preconditioner", "Ifpack"); 
  //lsParams.set("Preconditioner", "AztecOO"); 
  //lsParams.set("Aztec Preconditioner", "ilu"); 
  //lsParams.set("Overlap", 2);  
  //lsParams.set("Graph Fill", 2); 
  //lsParams.set("Aztec Preconditioner", "ilut"); 
  //lsParams.set("Overlap", 2);   
  //lsParams.set("Fill Factor", 2);   
  //lsParams.set("Drop Tolerance", 1.0e-12);   
  //lsParams.set("Aztec Preconditioner", "Polynomial"); 
  //lsParams.set("Polynomial Order", 6); 
#endif

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RCP<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));

  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RCP<Epetra_RowMatrix> A = 
    Teuchos::rcp(&Problem.getJacobian(),false);

#ifdef DO_XYZT
  Epetra_MultiVector initGuess(soln.Map(), globalComm->NumTimeStepsOnDomain());
  for (int i=0; i<globalComm->NumTimeStepsOnDomain(); i++) *(initGuess(i)) = soln;

#ifdef DO_XYZT_PREC
  // Sublist for linear solver of the preconditioner
  Teuchos::RCP<Teuchos::ParameterList> precLSParams =
       Teuchos::rcp(new Teuchos::ParameterList(lsParams));
  precLSParams->set("Aztec Solver", "GMRES");  
  precLSParams->set("Preconditioner", "Ifpack");  
  //precLSParams->set("Preconditioner", "AztecOO");  
  precLSParams->set("Max Iterations", 800);  
  precLSParams->set("Tolerance", 1e-6);
  precLSParams->set("Output Frequency", 50);    
  //precLSParams->set("XYZTPreconditioner", "None"); 
  precLSParams->set("XYZTPreconditioner", "Global"); 
  //precLSParams->set("XYZTPreconditioner", "Sequential"); 
  //precLSParams->set("XYZTPreconditioner", "Parallel"); 
  //precLSParams->set("XYZTPreconditioner", "Parareal"); 
  //precLSParams->set("XYZTPreconditioner", "BlockDiagonal"); 

  Teuchos::RCP<Teuchos::ParameterList> precPrintParams =
       Teuchos::rcp(new Teuchos::ParameterList(printParams));
  precPrintParams->set("MyPID", MyPID); 
  precPrintParams->set("Output Precision", 3);
  precPrintParams->set("Output Processor", 0);
  precPrintParams->set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  Teuchos::RCP<LOCA::Epetra::Interface::xyzt> ixyzt = 
    Teuchos::rcp(new LOCA::Epetra::Interface::xyzt(interface, 
						   initGuess, A,
						   globalComm, 
                                                   initCond, dt,
						   precPrintParams.get(), precLSParams.get()));

  Teuchos::RCP<Epetra_RowMatrix> Axyzt =
     Teuchos::rcp(&(ixyzt->getJacobian()),false);
  Epetra_Vector& solnxyzt = ixyzt->getSolution();
  Teuchos::RCP<Epetra_Operator> Mxyzt = 
     Teuchos::rcp(&(ixyzt->getPreconditioner()),false);

  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = ixyzt;

  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = ixyzt;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = 
    Teuchos::rcp(&(ixyzt->getPreconditioner()),false);
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iJac, Axyzt, iPrec, Mxyzt, 
						      solnxyzt));
#else
  Teuchos::RCP<LOCA::Epetra::Interface::xyzt> ixyzt = 
    Teuchos::rcp(new LOCA::Epetra::Interface::xyzt(interface, 
						   initGuess, A,
                                                   initCond, dt,
						   globalComm));

  Teuchos::RCP<Epetra_RowMatrix> Axyzt =
     Teuchos::rcp(&(ixyzt->getJacobian()),false);
  Epetra_Vector& solnxyzt = ixyzt->getSolution();

  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = ixyzt;

  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = ixyzt;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, Axyzt, 
						      solnxyzt));
#endif
  NOX::Epetra::Vector initialGuess(solnxyzt);
#else
  // Use an Epetra Scaling object if desired
  Teuchos::RCP<Epetra_Vector> scaleVec = 
    Teuchos::rcp(new Epetra_Vector(soln));
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);

  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, A, soln));
		                                      //&scaling);

  // Create the Group
  NOX::Epetra::Vector initialGuess(Teuchos::rcp(&soln,false), 
				   NOX::Epetra::Vector::CreateView,
				   NOX::DeepCopy);
#endif

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("alpha",0.6);
  pVector.addParameter("beta",2.0);

  // Create Epetra factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RCP<LOCA::GlobalData> globalData = 
    LOCA::createGlobalData(paramList, epetraFactory);

  Teuchos::RCP<LOCA::Epetra::Group> grp =
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams,
                 iReq, initialGuess, linSys, pVector));

  grp->computeF();

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					    NOX::StatusTest::NormF::Unscaled));
  //NOX::StatusTest::NormF relresid(*grp.get(), 1.0e-2);
  //NOX::StatusTest::NormUpdate update(1.0e-5);
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  //NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
  //converged.addStatusTest(absresid);
  //converged.addStatusTest(relresid);
  //converged.addStatusTest(wrms);
  //converged.addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(absresid);
  combo->addStatusTest(maxiters);

  // Create the method
  //NOX::Solver::Manager solver(grp, combo, nlParams);
  LOCA::Stepper stepper(globalData, grp, combo, paramList);

  // Initialize time integration parameters
#ifdef DO_XYZT
  int maxTimeSteps = 1; // No longer need a time integration loop
#else
  int maxTimeSteps = 2;
#endif
  int timeStep = 0;
  double time = 0.;
  
  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;
  
    globalData->locaUtils->out() 
      << "Time Step: " << timeStep << ",\tTime: " << time << std::endl;
  
//    NOX::StatusTest::StatusType status = solver.solve();
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status == LOCA::Abstract::Iterator::Finished)
      globalData->locaUtils->out() << "All tests passed" << std::endl;
    else
       globalData->locaUtils->out() << "Stepper failed to converge!" << std::endl;


    // Get the Epetra_Vector with the final solution from the solver
    const LOCA::Epetra::Group& finalGroup = 
      dynamic_cast<const LOCA::Epetra::Group&>(*(stepper.getSolutionGroup()));
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // End Nonlinear Solver **************************************

#ifndef DO_XYZT

    //Problem.reset(finalSolution);
    grp->setX(finalSolution);
    stepper.reset(globalData, grp, combo, paramList);
    grp->computeF();
#endif

  } // end time step while loop

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

  // Output timing info
  globalData->locaUtils->out() << "\nTimings :\n\tWallTime --> " << 
    myTimer.WallTime() - startWallTime << " sec."
	      << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
	      << " sec." << std::endl << std::endl;

  LOCA::destroyGlobalData(globalData);

return 0 ;
}

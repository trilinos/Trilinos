//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
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
#include "NOX.H"
#include "NOX_Epetra.H"

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

#ifdef HAVE_NOX_ML_EPETRA
#include "Teuchos_ParameterList.hpp"
#endif

// Headers needed for FD coloring 
#include <vector> 
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h" 
#endif

// New coupling library headers
#include "NOX_Multiphysics_Solver_Manager.H" 

// User's application specific files 
#include "Problem_Manager.H" 
#include "Problem_Interface.H" 
#include "Equation_A.H"              
#include "Equation_B.H"              
#include "Burgers.H"              

#include "HMX_PDE.H"              

#include "ConvDiff_PDE.H"              

// Added to allow timings
#include "Epetra_Time.h"

using namespace std;

int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Make sure we have a valid candidate command line
  if (argc < 2) 
  {
    cout << "Usage: " << argv[0] 
         << " -n number_of_elements [ -fixedPt -matlab -offBlocks -noFlow -hmx -convdiff ]" 
         << endl;
    exit(1);
  }

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Default run-time options that can be changed from the command line
  bool runMF         = true  ;
  bool useMatlab     = false ;
  bool doOffBlocks   = false ;
  bool noFlow        = true ;
  bool doBrusselator = true  ;  // default to Brusselator problem
  bool doHMX         = false ; 
  bool doConvDiff    = false ; 

  
  string inConcat("");

  // Dump out the command line args
  for( int i = 1; i < argc; ++i )
  {
    cout << "arg [" << i << "] = " << argv[i] << endl;
    inConcat += argv[i];
    inConcat += " ";
  }
  
  inConcat.erase( inConcat.size(), 1 );

  cout << "ECHO = " << inConcat << endl;

  string option = "-fixedPt";
  if( inConcat.find(option) != string::npos )
  {
    runMF = false;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-matlab";
  if( inConcat.find(option) != string::npos )
  {
    useMatlab = true;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-offBlocks";
  if( inConcat.find(option) != string::npos )
  {
    doOffBlocks = true;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-noFlow";
  if( inConcat.find(option) != string::npos )
  {
    noFlow = true;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-hmx";
  if( inConcat.find(option) != string::npos )
  {
    doBrusselator = false;
    doHMX         = true ;
    doConvDiff    = false;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-convdiff";
  if( inConcat.find(option) != string::npos )
  {
    doBrusselator = false;
    doHMX         = false;
    doConvDiff    = true ;
    inConcat.replace( inConcat.find(option), option.size()+1, "");
  }

  option = "-n";
  if( inConcat.find(option) != string::npos )
  {
    inConcat.replace( inConcat.find(option), option.size()+1, "");
    //cout << "Found number of elements.\n" << "Remaining command line :\n " << inConcat << endl;
  }

  int NumGlobalNodes = atoi(inConcat.c_str()) + 1;

  // The number of unknowns must be at least equal to the number of processors.
  if (NumGlobalNodes < NumProc) 
  {
    cout << "numGlobalNodes = " << NumGlobalNodes 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Begin Nonlinear Solver ************************************

  // NOTE: For now these parameters apply to all problems handled by
  // Problem_Manager.  Each problem could be made to have its own
  // parameter list as wwell as its own convergence test(s).

  // Create the top level parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

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
			NOX::Utils::Warning);

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
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 50);    
  //lsParams.set("Preconditioning", "None");   
  //lsParams.set("Preconditioning", "AztecOO: Jacobian Matrix");   
  lsParams.set("Preconditioner", "AztecOO");   
  //lsParams.set("Preconditioner", "Ifpack");
  //lsParams.set("Graph Fill", 2);
  //lsParams.set("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.set("Preconditioning", "User Supplied Preconditioner");
  //lsParams.set("Aztec Preconditioner", "ilu"); 
  //lsParams.set("Overlap", 2);  
  //lsParams.set("Graph Fill", 2); 
  //lsParams.set("Aztec Preconditioner", "ilut"); 
  //lsParams.set("Overlap", 2);   
  //lsParams.set("Fill Factor", 2);   
  //lsParams.set("Drop Tolerance", 1.0e-12);   
  //lsParams.set("Aztec Preconditioner", "Polynomial"); 
  //lsParams.set("Polynomial Order", 6); 
#ifdef HAVE_NOX_ML_EPETRA
  
  //lsParams.set("Preconditioner", "ML");
  //Teuchos::ParameterList MLList;
  //if( lsParams.get("Preconditioner", "None") == "ML" ) {
  //  // This Teuchos parameter list is needed for ML
  //
  //  // These specifications come straight from the example in 
  //  // Trilinos/packages/ml/example/ml_example_epetra_preconditioner.cpp
  //
  //  // set defaults for classic smoothed aggregation
  //  ML_Epetra::SetDefaults("SA",MLList);
  //  // maximum number of levels
  //  MLList.set("max levels",5);
  //  MLList.set("increasing or decreasing","decreasing");
  //  // use Uncoupled scheme to create the aggregate,
  //  // from level 3 use the better but more expensive MIS
  //  MLList.set("aggregation: type", "Uncoupled");
  //  MLList.set("aggregation: type (level 3)", "MIS");
  //  // smoother is Gauss-Seidel. Example file
  //  // ml_example_epetra_preconditioner_2level.cpp shows how to use
  //  // AZTEC's preconditioners as smoothers
  //  MLList.set("smoother: type","Gauss-Seidel");
  //  // use both pre and post smoothing
  //  MLList.set("smoother: pre or post", "both");
  //  // solve with serial direct solver KLU
  //  MLList.set("coarse: type","Jacobi");
  //
  //  lsParams.set("ML Teuchos Parameter List", &MLList);
  //}
#endif

  // Create the convergence tests
  // Note: as for the parameter list, both (all) problems use the same 
  // convergence test(s) for now, but each could have its own.
  Teuchos::RefCountPtr<NOX::StatusTest::NormF>          absresid  =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate>     update    = 
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo>          converged = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  //converged->addStatusTest(update);
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> finiteValue = 
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteValue);

  // Create the Problem Manager
  Problem_Manager problemManager(Comm, doOffBlocks, 0, useMatlab);

  // Note that each problem could contain its own nlParams list as well as
  // its own convergence test(s). 
  problemManager.registerParameters(nlParamsPtr);
  problemManager.registerStatusTest(combo);

  // Allow one of two supported tests
  if( doBrusselator ) 
  {
    // Create each part of the Brusselator problem class.  
    Equation_A ProblemA(Comm, NumGlobalNodes, "Temperature" );
    Equation_B ProblemB(Comm, NumGlobalNodes, "Species"     );
    Burgers  burgers   (Comm, NumGlobalNodes, "Burgers"     );
  
    // An interesting note: the order of solving each problem is based on the
    // order of adding.  For this decoupled problem, problem B is linear
    // with respect to its variables, whereas problem A is nonlinear wrt to its
    // variables.  The order of solution appears to strongly affect the rate
    // of convergence of the decoupled Brusselator.  Solving problem A first
    // dramatically reduces the number of total iterations.
    problemManager.addProblem(ProblemA);
    problemManager.addProblem(ProblemB);
    if( !noFlow )
      problemManager.addProblem(burgers);
  
  //  problemManager.createDependency("Temperature", "Species");
    problemManager.createDependency(ProblemA, ProblemB);
    problemManager.createDependency(ProblemB, ProblemA);

    if( !noFlow )
    {
      problemManager.createDependency(ProblemA, burgers);
      problemManager.createDependency(ProblemB, burgers);
      problemManager.createDependency(burgers, ProblemA);
    }
  
    problemManager.registerComplete(); // Trigger setup of groups, solvers, etc.
  
    problemManager.outputStatus(std::cout);
  
    // Initialize time integration parameters
    int maxTimeSteps = 3;
    int timeStep = 0;
    double time = 0.;
    double dt = 0.100;
    problemManager.setAlldt(dt);
    
    // Print initial solution
    char file_name[25];
    FILE *ifp;
    Epetra_Vector& xMesh = ProblemA.getMesh();
    int NumMyNodes = xMesh.Map().NumMyElements();
    (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (int i=0; i<NumMyNodes; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
                                   xMesh[i], (*ProblemA.getSolution())[i], 
                                   (*ProblemB.getSolution())[i]);
    fclose(ifp);
    //FILE *ifp2;
    //Epetra_Vector& burgersX = burgers.getMesh();
    //(void) sprintf(file_name, "burgers.%d_%d",MyPID,timeStep);
    //ifp2 = fopen(file_name, "w");
    //for (int i = 0; i < 10*NumMyNodes; ++i)
    //  fprintf(ifp2, "%d  %E  %E\n", burgersX.Map().MinMyGID()+i, 
    //                               burgersX[i], burgers.getSolution()[i]);
    //fclose(ifp2);
    
    
    // Time integration loop
    while(timeStep < maxTimeSteps) {
  
      timeStep++;
      time += dt;
    
      cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
    
      // Solve the coupled problem
      if( runMF )
        problemManager.solveMF(); // Need a status test check here ....
      else
      {
        // Create the loose coupling solver manager
        Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> > solversVec =
          Teuchos::rcp( new vector<NOX::Solver::Manager*> );

        map<int, NOX::Solver::Manager*>::iterator iter = problemManager.getSolvers().begin(),
                                              iter_end = problemManager.getSolvers().end()   ;
        for( ; iter_end != iter; ++iter )
        {
          cout << " ........  registered Solver::Manager # " << (*iter).first << endl;
          solversVec->push_back( (*iter).second );
        }

        // Package the Problem_Manager as the DataExchange::Intreface
        Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface> dataExInterface =
          Teuchos::rcp( &problemManager, false );

        NOX::Multiphysics::Solver::Manager cplSolv( solversVec, dataExInterface, combo, nlParamsPtr );

        cplSolv.solve();
      }

      problemManager.outputSolutions(timeStep);
  
      // Reset problems by copying solution into old solution
      problemManager.resetProblems();
  
    } // end time step while loop
  
    // Output timing info
    if(MyPID==0)
      cout << "\nTimings :\n\tWallTime --> " << 
  	    myTimer.WallTime() - startWallTime << " sec."
           << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
           << " sec." << endl << endl;

  }
  else if( doHMX )
  { // HMX Cook-off problem

    string nameT 		= "Temperature";
    double Const_R		= 1.9872 ;
    double Specific_H		= 0.42 ;
    double Density		= 1.90 ;
    double Thermal_K		= 0.8658e-3 ;
    double diffCoef_T		= Thermal_K / (Density * Specific_H);

    double StericCoef_T		= 0.0;
    double PreExp_T 		= 0.0;
    double ActEnergy_T	= 0.0 ;
    map<string, double> SrcTermExponent_T; // Leave empty if no volume source
    map<string, double> SrcTermWeight_T; 
      SrcTermWeight_T.insert( pair<string, double> ("SpeciesA", -190.0) );
      SrcTermWeight_T.insert( pair<string, double> ("SpeciesB",  570.0) );
      SrcTermWeight_T.insert( pair<string, double> ("SpeciesC", 2280.0) );

    // Create each part of the HMX cook-off problem
    HMX_PDE HMX_TempEq (Comm, 
                    diffCoef_T,
                    Const_R,
		    StericCoef_T, // Dummy for Temp Eq.
		    PreExp_T, // Dummy for Temp Eq.
		    ActEnergy_T, // Dummy for Temp Eq.
		    SrcTermExponent_T,
		    SrcTermWeight_T,
		    1*NumGlobalNodes, nameT);

    HMX_TempEq.setTempFieldName(HMX_TempEq.getName());

    // Override default initialization with values we want
    HMX_TempEq.initializeSolution(500.0);

    problemManager.addProblem(HMX_TempEq);


    string nameA 		= "SpeciesA";
    double diffCoef_A 		= 0.0 ;
    //double stericCoef_A		= 0.0 ;
    double stericCoef_A		= 2.0 ; // ROGER: high coupling
    //double stericCoef_A		= 2.0 ; // ROGER: high coupling
    double preExp_A 		= exp(48.7) ;
    double actEnergy_A 		= 52700.0 ;
    map<string, double> SrcTermExponent_A; 
      SrcTermExponent_A.insert( pair<string, double> (nameA, 1.0) );
    map<string, double> SrcTermWeight_A; 
      SrcTermWeight_A.insert( pair<string, double> (nameA, -1.0) );

    HMX_PDE HMX_RxnA(Comm, 
                  diffCoef_A,
                  Const_R,
                  stericCoef_A,
                  preExp_A,
                  actEnergy_A,
                  SrcTermExponent_A,
                  SrcTermWeight_A,
                  1*NumGlobalNodes, nameA);

    HMX_RxnA.setTempFieldName(HMX_TempEq.getName());

    // Override default initialization with values we want
    HMX_RxnA.initializeSolution(2.0);

    problemManager.addProblem(HMX_RxnA);


    string nameB 		= "SpeciesB";
    double diffCoef_B 		= 0.0 ;
    double stericCoef_B		= 1.0 ;
    //double stericCoef_B		= -1.0 ;
    double preExp_B 		= exp(37.5) ;
    double actEnergy_B 		= 44100.0 ;
    map<string, double> SrcTermExponent_B; 
      SrcTermExponent_B.insert( pair<string, double> (nameB, 1.0) );
    map<string, double> SrcTermWeight_B; 
      SrcTermWeight_B.insert( pair<string, double> (nameA, 1.0) );
      SrcTermWeight_B.insert( pair<string, double> (nameB, -1.0) );

    HMX_PDE HMX_RxnB(Comm, 
                  diffCoef_B,
                  Const_R,
                  stericCoef_B,
                  preExp_B,
                  actEnergy_B,
                  SrcTermExponent_B,
                  SrcTermWeight_B,
                  NumGlobalNodes, nameB);

    HMX_RxnB.setTempFieldName(HMX_TempEq.getName());

    // Override default initialization with values we want
    HMX_RxnB.initializeSolution(1.0);

    problemManager.addProblem(HMX_RxnB);
  

    string nameC 		= "SpeciesC";
    double diffCoef_C 		= 0.0 ;
    double stericCoef_C		= 0.0 ;
    double preExp_C 		= exp(28.1) ;
    double actEnergy_C 		= 34100.0 ;
    map<string, double> SrcTermExponent_C; 
      SrcTermExponent_C.insert( pair<string, double> (nameC, 2.0) );
    map<string, double> SrcTermWeight_C; 
      SrcTermWeight_C.insert( pair<string, double> (nameB, 2.0) );
      SrcTermWeight_C.insert( pair<string, double> (nameC, -2.0) );

    HMX_PDE HMX_RxnC(Comm, 
                  diffCoef_C,
                  Const_R,
                  stericCoef_C,
                  preExp_C,
                  actEnergy_C,
                  SrcTermExponent_C,
                  SrcTermWeight_C,
                  1*NumGlobalNodes, nameC);

    HMX_RxnC.setTempFieldName(HMX_TempEq.getName());

    // Override default initialization with values we want
//    HMX_RxnC.initializeSolution(0.0);

    problemManager.addProblem(HMX_RxnC);
  
    problemManager.createDependency(HMX_TempEq, HMX_RxnA);
    problemManager.createDependency(HMX_TempEq, HMX_RxnB);
    problemManager.createDependency(HMX_TempEq, HMX_RxnC);

    problemManager.createDependency(HMX_RxnA, HMX_TempEq);

    problemManager.createDependency(HMX_RxnB, HMX_TempEq);
    problemManager.createDependency(HMX_RxnB, HMX_RxnA);

    problemManager.createDependency(HMX_RxnC, HMX_TempEq);
    problemManager.createDependency(HMX_RxnC, HMX_RxnB);
  
    problemManager.registerComplete();
  
    problemManager.outputStatus(std::cout);
  
    // Initialize time integration parameters
    int maxTimeSteps = 5;
    int timeStep = 0;
    double time = 0.;
    //double dt = HMX_TempEq.getdt();
    double dt = 10.0 * HMX_TempEq.getdt();
    
    // Print initial solution
    char file_name[25];
    FILE *ifp;
    Epetra_Vector& xMesh = HMX_TempEq.getMesh();
    int NumMyNodes = xMesh.Map().NumMyElements();
    (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (int i=0; i<NumMyNodes; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, 
                                   xMesh[i], (*HMX_TempEq.getSolution())[i], 
                                   (*HMX_RxnA.getSolution())[i]);
    fclose(ifp);
    
    // Time integration loop
    while(timeStep < maxTimeSteps) {
  
      timeStep++;
      time += dt;
    
      cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
    
      // Solve the coupled problem
      if( runMF )
        problemManager.solveMF(); // Need a status test check here ....
      else
        problemManager.solve(); // Need a status test check here ....
    
      problemManager.outputSolutions(timeStep);
  
      // Reset problems by copying solution into old solution
      problemManager.resetProblems();
  
    } // end time step while loop
  
    // Output timing info
    if(MyPID==0)
      cout << "\nTimings :\n\tWallTime --> " << 
  	    myTimer.WallTime() - startWallTime << " sec."
           << "\n\tElapsedTime --> " << myTimer.ElapsedTime() 
           << " sec." << endl << endl;
  }
  else if( doConvDiff )
  { // Convection-Diffusion coupled problem

    // Coupling parameters
    double alpha		= 0.50          ;
    double beta 		= 0.40          ;

    // Domain boundary temperaures
    double Tleft  		= 0.98          ;
    double Tright 		= 1.0           ;

    // Distinguish certain parameters needed for T1_analytic
    double peclet_1     	= 9.0           ;
    double peclet_2     	= 0.0           ;
    double kappa_1      	= 1.0           ;
    double kappa_2		= 0.1           ;

    //double T1_analytic          = ( kappa_2*(1.0-exp(peclet_1))*Tright - Tleft*peclet_1*exp(peclet_1) ) /
    //                              ( kappa_2*(1.0-exp(peclet_1)) - peclet_1*exp(peclet_1) );
    //double T1_analytic          = 0.99430092; // 5.67
    double T1_analytic          = 0.99050495; // 2.50
    //double T1_analytic          = 0.98247093; // 0.3
    //double T1_analytic          = 0.98177822; // 0.2

    // Create Region 1 PDE
    string myName 		= "Region_1"    ;
    double radiation		= 0.0           ;
    double xmin  		= 0.0           ;
    double xmax  		= 1.0           ;

    ConvDiff_PDE Reg1_PDE (
                    Comm, 
                    peclet_1,
                    radiation,
		    kappa_1,
		    alpha,
		    xmin,
		    xmax, 
                    Tleft,
                    T1_analytic,
		    1*NumGlobalNodes, 
                    myName  );

    // Override default initialization with values we want
    Reg1_PDE.initializeSolution(0.995);

    problemManager.addProblem(Reg1_PDE);


    // Create Region 2 PDE
    myName 		        = "Region_2"    ;
    radiation		        = 5.67          ;
    xmin  		        = 1.0           ;
    xmax  		        = 2.0           ;

    ConvDiff_PDE Reg2_PDE (
                    Comm, 
                    peclet_2,
                    radiation,
		    kappa_2,
		    beta,
		    xmin,
		    xmax, 
                    T1_analytic,
                    Tright,
		    1*NumGlobalNodes, 
                    myName  );

    // Override default initialization with values we want
    Reg2_PDE.initializeSolution(0.995);

    problemManager.addProblem(Reg2_PDE);

    problemManager.createDependency(Reg1_PDE, Reg2_PDE);
    problemManager.createDependency(Reg2_PDE, Reg1_PDE);
  
    problemManager.registerComplete();
  
    problemManager.outputStatus(std::cout);

    cout << "\n\tAnalytic solution, T_1 = " << T1_analytic << "\n" << endl;
  
    // Print initial solution
    problemManager.outputSolutions(0);

    // Solve the coupled problem
    if( runMF )
      problemManager.solveMF(); // Need a status test check here ....
    else
      problemManager.solve(); // Need a status test check here ....
    
    problemManager.outputSolutions(1);
  
    // Output timing info
    if(MyPID==0)
      cout << "\nTimings :\n\tWallTime --> " << myTimer.WallTime() - startWallTime << " sec."
           << "\n\tElapsedTime --> " << myTimer.ElapsedTime() << " sec." << endl << endl;
  }
  else
  {
    if(MyPID==0)
      cout << "Test failed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

    return -1;
  }

  // Need to put in a check for convergence
  // Added the following so test actually passes in parallel
  if(MyPID==0)
    cout << "Test passed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
}

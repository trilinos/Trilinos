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
                                                                    
// 1D Transient Finite Element Test Problem
/* Solves the nonlinear PDE for u(x,t):
*
* du     d2u    8
* --- =  --- + --- u**2( 1 - u )
* dt     dx2   k*k
*
* subject to u(-infinity,t) = 1
*            u(-infinity,t) = 0
* and
*            u(x,0) = 0.5 * ( 1 - tanh(x/k) )
*
* with d representing partial differentiation.
*
* An exact closed solution is the following:
*
*                             x - 2t/k
*  u(x,t) = 0.5 * ( 1 - tanh(--------- ) )
*                                k
*
* This problem is examined with a variety of time integration schemes in:
* "Studies on the Convergence of Various Time-Integration Schemes for the
* Radiation-Diffusion Problem," Curtis C. Ober & John N. Shadid, in prep.
*
* In this example, only a 1st-order fully implicit (backward Euler)
* time integration scheme is considered currently.
*
* Values for k, time step size, and finite spatial extent are specified in
* the constructor initialization list in FiniteElementProblem.C using
* variables factor, dt and xmin,xmax, respectively.
* The number of time steps to be taken is specified by variable
* maxTimeSteps below.
*/

#undef PRINT_RESULTS_TO_FILES

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

// User's application specific files 
#include "Interface_EpetraNew.H" 

using namespace std;

// ------------------------------------------------------------------------
// ------------------  Declaration with definition ------------------------
// ------------------------------------------------------------------------
class TransientInterface : public Interface 
{

public:

  // ---------------------------
  // --------------- Constructor ------
  // ---------------------------
  TransientInterface(int NumGlobalElements, Epetra_Comm& Comm, 
                     double xmin = 0.0, double xmax = 1.0) :
    Interface(NumGlobalElements, Comm, xmin, xmax),
    oldSolution(0),
    exactSolution(0)
    {
      // Must reinitialize the solution now that the inherited class exists
      initializeSoln();

      oldSolution = new Epetra_Vector(*initialSolution);
    };

  // ---------------------------
  // --------------- Destructor ------
  // ---------------------------
  virtual ~TransientInterface()
    { 
      delete oldSolution;
      if (exactSolution) delete exactSolution;
    };

  // ---------------------------
  // --------------- Matrix and Residual Fills
  // ---------------------------
  bool evaluate(NOX::EpetraNew::Interface::Required::FillType flag, 
                const Epetra_Vector* soln, 
  		Epetra_Vector* tmp_rhs, 
  		Epetra_RowMatrix* tmp_matrix)
  {

    //Determine what to fill (F or Jacobian)
    bool fillF = false;
    bool fillMatrix = false;
    if (tmp_rhs != 0) {
      fillF = true;
      rhs = tmp_rhs;
    }
    else {
      fillMatrix = true;
    }
  
    // "flag" can be used to determine how accurate your fill of F should be
    // depending on why we are calling evaluate (Could be using computeF to 
    // populate a Jacobian or Preconditioner).
    if (flag == NOX::EpetraNew::Interface::Required::Residual) {
      // Do nothing for now
    }
    else if (flag == NOX::EpetraNew::Interface::Required::Jac) {
      // Do nothing for now
    }
    else if (flag == NOX::EpetraNew::Interface::Required::Prec) {
      // Do nothing for now
    }
    else if (flag == NOX::EpetraNew::Interface::Required::User) {
      // Do nothing for now
    }
  
  
    // Create the overlapped solution and position vectors
    Epetra_Vector u(*OverlapMap);
    Epetra_Vector uold(*OverlapMap);
    Epetra_Vector xvec(*OverlapMap);
  
    // Export Solution to Overlap vector
    u.Import(*soln, *Importer, Insert);
    uold.Import(*oldSolution, *Importer, Insert);
    xvec.Import(*xptr, *Importer, Insert);
  
    // Declare required variables
    int ierr;
    int OverlapNumMyElements = OverlapMap->NumMyElements();
  
    int OverlapMinMyGID;
    if (MyPID == 0) OverlapMinMyGID = StandardMap->MinMyGID();
    else OverlapMinMyGID = StandardMap->MinMyGID()-1;
  
    int row, column;
    double jac;
    double xx[2];
    double uu[2];
    double uuold[2];
    Basis basis;
  
    // Zero out the objects that will be filled
    if (fillF) 
      rhs->PutScalar(0.0);
    if (fillMatrix) 
      Jacobian->PutScalar(0.0);
  
    // Loop Over # of Finite Elements on Processor
    for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
      
      // Loop Over Gauss Points
      for(int gp=0; gp < 2; gp++) {
        // Get the solution and coordinates at the nodes 
        xx[0]=xvec[ne];
        xx[1]=xvec[ne+1];
        uu[0]=u[ne];
        uu[1]=u[ne+1];
        uuold[0]=uold[ne];
        uuold[1]=uold[ne+1];
        // Calculate the basis function at the gauss point
        basis.computeBasis(gp, xx, uu, uuold);
  	            
        // Loop over Nodes in Element
        for (int i=0; i< 2; i++) {
  	row=OverlapMap->GID(ne+i);
  	//printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
  	//     MyPID, row, ne+i,StandardMap.MyGID(row));
  	if (StandardMap->MyGID(row)) {
  	  if (fillF) {
            (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
              +basis.wt*basis.dx*(
                (basis.uu-basis.uuold)/dt*basis.phi[i] +
                (1.0/(basis.dx*basis.dx))*basis.duu*basis.dphide[i] -
                8.0/factor/factor*basis.uu*basis.uu*
                (1.0-basis.uu)*basis.phi[i]);
  	  }
  	}
  	// Loop over Trial Functions
  	if (fillMatrix) {
  	  for(int j=0;j < 2; j++) {
  	    if (StandardMap->MyGID(row)) {
  	      column=OverlapMap->GID(ne+j);
              jac=basis.wt*basis.dx*(
                basis.phi[j]/dt*basis.phi[i] +
	        (1.0/(basis.dx*basis.dx))*
                basis.dphide[j]*basis.dphide[i] -
                8.0/factor/factor*
                (2.0*basis.uu-3.0*basis.uu*basis.uu)*
                basis.phi[j]*basis.phi[i]);
  	      ierr=Jacobian->SumIntoGlobalValues(row, 1, &jac, &column);
  	    }
  	  }
  	}
        }
      }
    } 
  
    // Insert Boundary Conditions and modify Jacobian and function (F)
    // U(0)=1
    if (MyPID==0) {
      if (fillF) 
        (*rhs)[0]= (*soln)[0] - 1.0;
      if (fillMatrix) {
        column=0;
        jac=1.0;
        Jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
        column=1;
        jac=0.0;
        Jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
      }
    }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(xmax)=0
  if (MyPID==NumProc-1) {
    if (fillF)
      (*rhs)[NumMyElements-1]= (*soln)[OverlapNumMyElements-1] - 0.0;
    if (fillMatrix) {
      row=NumGlobalElements-1;
      column=row;
      jac=1.0;
      Jacobian->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      jac=0.0;
      Jacobian->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

    // Sync up processors to be safe
    Comm->Barrier();
   
    Jacobian->TransformToLocal();
  
    return true;
  }

  // ---------------------------
  // --------------- Set desired initial condition / initial guess
  // ---------------------------
  bool initializeSoln()
  {
    Epetra_Vector& x = *xptr;
  
    double arg;
    for(int i=0; i<NumMyElements; i++) {
      arg = x[i]/factor;
        (*initialSolution)[i] = (1.0 - ( exp(arg) - exp(-arg) ) /
        ( exp(arg) + exp(-arg) )) / 2.0;
    }
    return true;
  }

  // Reset problem for next parameter (time) step.
  // For now, this simply updates oldsoln with the given Epetra_Vector
  bool reset(const Epetra_Vector& x) 
    { *oldSolution = x;  return true; };

  // Return a reference to the Epetra_Vector with the old solution
  Epetra_Vector& getOldSoln()
    { return *oldSolution; };

  // ---------------------------
  // --------------- Return a reference to the Epetra_Vector 
  // --------------- with the exact analytic solution
  // ---------------------------
  Epetra_Vector& getExactSoln(double time)
    {  
      // Create the exact solution vector only if needed and if not already
      // created.
      if( !exactSolution )
        exactSolution = new Epetra_Vector(*xptr);

      Epetra_Vector& x = *xptr;

      for(int i=0; i<NumMyElements; i++)
        (*exactSolution)[i] = 
          (1.0 - tanh( (x[i]-2.0*time/factor)/factor )) / 2.0;

      return *exactSolution;
    }

  // Accesor function for setting time step
  void setdt( double dt_ ) { dt = dt_; }
  
  // Accesor function for time obtaining step
  double getdt() { return dt; }
  
private:

  double dt; 	// Time step size

  Epetra_Vector *oldSolution;
  Epetra_Vector *exactSolution;

};


// ------------------------------------------------------------------------
// ---------------------------   Main Program -----------------------------
// ------------------------------------------------------------------------
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

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  bool verbose = false;
  // Check for verbose output
  if (argc>1) 
    if (argv[1][0]=='-' && argv[1][1]=='v') 
      verbose = true;

  // Get the number of elements from the command line
  int NumGlobalElements = 0;
  if ((argc > 2) && (verbose))
    NumGlobalElements = atoi(argv[2]) + 1;
  else if ((argc > 1) && (!verbose))
    NumGlobalElements = atoi(argv[1]) + 1;
  else 
    NumGlobalElements = 101;

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    throw "NOX Error";
  }

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  TransientInterface interface(NumGlobalElements, Comm, -20.0, 20.0);

  // Set the PDE nonlinear coefficient for this problem
  interface.setPDEfactor(1.0);

  // Get the vector from the Problem
  Epetra_Vector& soln = interface.getSolution();

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  if (verbose)
    printParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning +
                             NOX::Utils::Debug +
			     NOX::Utils::Error);
  else
    printParams.setParameter("Output Information", NOX::Utils::Error);

  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);  
  //lsParams.setParameter("Preconditioner", "None");
  lsParams.setParameter("Preconditioner", "AztecOO");
  //lsParams.setParameter("Jacobian Operator", "Finite Difference");
  //lsParams.setParameter("Jacobian Operator", "Matrix-Free");
  //lsParams.setParameter("Preconditioner Operator", "Finite Difference");

  // Use an Epetra Scaling object if desired
  //Epetra_Vector scaleVec(soln);
  //NOX::Epetra::Scaling scaling;
  //scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  // Create all possible Epetra_Operators.
  // 1. User supplied (Epetra_RowMatrix)
  Epetra_RowMatrix& Analytic = interface.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  NOX::EpetraNew::MatrixFree MF(interface, soln);
  // 3. Finite Difference (Epetra_RowMatrix)
  NOX::EpetraNew::FiniteDifference FD(interface, soln);

  // Four constructors to create the Linear System
  NOX::EpetraNew::Interface::Required& iReq = interface;

  // **** Ctor #1 - No Jac and No Prec
  //NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams, 
  //				      iReq, soln);

  // **** Ctor #2 - Jac but no Prec
  NOX::EpetraNew::Interface::Jacobian& iJac = interface;
  NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
  				     iReq, iJac, Analytic, soln);

  // **** Ctor #3 - Prec but no Jac
  //NOX::EpetraNew::Interface::Preconditioner& iPrec = FD;
  //NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
  //				      iReq, iPrec, FD, soln);

  // **** Ctor #4 - Prec and Jac
  //NOX::EpetraNew::Interface::Jacobian& iJac = interface;
  //NOX::EpetraNew::Interface::Preconditioner& iPrec = interface;
  //NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
//					     iJac, Analytic, iPrec, Analytic, soln);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
  NOX::EpetraNew::Group grp(printParams, iReq, initialGuess, linSys);  

  // uncomment the following for loca supergroups
  MF.setGroupForComputeF(grp);
  FD.setGroupForComputeF(grp);

  // Test group accessor
  //const NOX::EpetraNew::LinearSystem& lins = grp.getLinearSystem();
  //lins.createPreconditioner(soln, lsParams, false);

  // Create the convergence tests
  NOX::StatusTest::NormF absresid(1.0e-8);
  NOX::StatusTest::NormF relresid(grp, 1.0e-2);
  NOX::StatusTest::NormUpdate update(1.0e-5);
  NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
  converged.addStatusTest(absresid);
  converged.addStatusTest(relresid);
  converged.addStatusTest(wrms);
  converged.addStatusTest(update);
  NOX::StatusTest::MaxIters maxiters(20);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  NOX::StatusTest::FiniteValue fv;
  combo.addStatusTest(fv);
  combo.addStatusTest(converged);
  combo.addStatusTest(maxiters);

  // Initialize time integration parameters
  int maxTimeSteps = 10;
  int timeStep = 0;
  double time = 0.0;
  double dt = 0.10;
  interface.setdt(dt);

#ifdef PRINT_RESULTS_TO_FILES
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E  %E\n", soln.Map().MinMyGID()+i,
             interface.getMesh()[i], soln[i]);
  fclose(ifp);
#endif

  // Create the solver
  NOX::Solver::Manager solver(grp, combo, nlParams);

  // Create a print class for controlling output below
  NOX::Utils utils(printParams);

  // Overall status flag
  int ierr = 0;

  // Time integration loop
  while(timeStep < maxTimeSteps) {
    timeStep++;
    time += dt;

    cout << "Time Step: " << timeStep << ",\tTime: " << time << endl;
  
    NOX::StatusTest::StatusType status = solver.solve();

    // Check for convergence
    if (status != NOX::StatusTest::Converged) {
        ierr++;
        if (utils.isPrintProcessAndType(NOX::Utils::Error))
          cout << "Nonlinear solver failed to converge!" << endl;
    }

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::EpetraNew::Group& finalGroup = dynamic_cast<const NOX::EpetraNew::Group&>(solver.getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    Epetra_Vector& exactSolution = interface.getExactSoln(time);

    // End Nonlinear Solver **************************************

#ifdef PRINT_RESULTS_TO_FILES
    // Print solution
    (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (int i=0; i<NumMyElements; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", soln.Map().MinMyGID()+i,
              interface.getMesh()[i], finalSolution[i],exactSolution[i]);
    fclose(ifp);
#endif

    interface.reset(finalSolution);
    grp.setX(finalSolution);
    solver.reset(grp, combo, nlParams);
    grp.computeF();

  } // end time step while loop

  // Output the parameter list
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
    cout << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ierr;
}



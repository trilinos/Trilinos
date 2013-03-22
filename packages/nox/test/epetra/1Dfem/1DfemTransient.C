//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
#include "1DfemInterface.H" 

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
  bool evaluate(NOX::Epetra::Interface::Required::FillType flag, 
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
    if (flag == NOX::Epetra::Interface::Required::Residual) {
      // Do nothing for now
    }
    else if (flag == NOX::Epetra::Interface::Required::Jac) {
      // Do nothing for now
    }
    else if (flag == NOX::Epetra::Interface::Required::Prec) {
      // Do nothing for now
    }
    else if (flag == NOX::Epetra::Interface::Required::User) {
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
      jacobian->PutScalar(0.0);
  
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
  	      ierr=jacobian->SumIntoGlobalValues(row, 1, &jac, &column);
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
        jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
        column=1;
        jac=0.0;
        jacobian->ReplaceGlobalValues(0, 1, &jac, &column);
      }
    }

  // Insert Boundary Conditions and modify Jacobian and function (F)
  // U(xmax)=0
  if (MyPID==NumProc-1) {
    if (fillF)
      (*rhs)[NumMyElements-1]= (*soln)[NumMyElements-1] - 0.0;
    if (fillMatrix) {
      row=NumGlobalElements-1;
      column=row;
      jac=1.0;
      jacobian->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      jac=0.0;
      jacobian->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

    // Sync up processors to be safe
    Comm->Barrier();
   
    jacobian->FillComplete();
  
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
    std::cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << std::endl;
    throw "NOX Error";
  }

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Teuchos::RCP<TransientInterface> interface = 
    Teuchos::rcp(new TransientInterface(NumGlobalElements, Comm, -20.0, 20.0));
  double dt = 0.10;
  interface->setdt(dt);

  // Set the PDE nonlinear coefficient for this problem
  interface->setPDEfactor(1.0);

  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = interface->getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  if (verbose)
    printParams.set("Output Information", 
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
    printParams.set("Output Information", NOX::Utils::Error);

  // Create a print class for controlling output below
  NOX::Utils utils(printParams);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 800);  
  lsParams.set("Tolerance", 1e-4);  
  lsParams.set("Preconditioner", "AztecOO");

  // Create all possible Epetra_Operators.
  // 1. User supplied (Epetra_RowMatrix)
  Teuchos::RCP<Epetra_RowMatrix> Analytic = interface->getJacobian();
  // 2. Matrix-Free (Epetra_Operator)

  // Four constructors to create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, Analytic, 
						      noxSoln));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					noxSoln, 
					linSys));  
  NOX::Epetra::Group& grp = *(grpPtr.get());

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Initialize time integration parameters
  int maxTimeSteps = 10;
  int timeStep = 0;
  double time = 0.0;

#ifdef PRINT_RESULTS_TO_FILES
  // Print initial solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E  %E\n", soln.Map().MinMyGID()+i,
             interface->getMesh()[i], soln[i]);
  fclose(ifp);
#endif

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

  // Overall status flag
  int ierr = 0;

  // Time integration loop
  while(timeStep < maxTimeSteps) {
    timeStep++;
    time += dt;

    utils.out() << "Time Step: " << timeStep << ",\tTime: " << time << std::endl;
  
    NOX::StatusTest::StatusType status = solver->solve();

    // Check for convergence
    if (status != NOX::StatusTest::Converged) {
        ierr++;
        if (utils.isPrintType(NOX::Utils::Error))
          utils.out() << "Nonlinear solver failed to converge!" << std::endl;
    }

    
    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    //Epetra_Vector& exactSolution = interface->getExactSoln(time);
    

    // End Nonlinear Solver **************************************

#ifdef PRINT_RESULTS_TO_FILES
    // Print solution
    (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
    ifp = fopen(file_name, "w");
    for (int i=0; i<NumMyElements; i++)
      fprintf(ifp, "%d  %E  %E  %E\n", soln.Map().MinMyGID()+i,
              interface->getMesh()[i], finalSolution[i],exactSolution[i]);
    fclose(ifp);
#endif

    interface->reset(finalSolution);
    grp.setX(finalSolution);
    solver->reset(grp.getX(), combo);
    grp.computeF();

  } // end time step while loop

  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << std::endl << "Final Parameters" << std::endl
	 << "****************" << std::endl;
    solver->getList().print(utils.out());
    utils.out() << std::endl;
  }

  // Test for convergence
  
#ifndef HAVE_MPI 
  // 1. Linear solve iterations on final time step (30)- SERIAL TEST ONLY!
  //    The number of linear iterations changes with # of procs.
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations",0) != 30) {
    ierr = 1;
  }
#endif
  // 2. Nonlinear solve iterations on final time step (3)
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) != 3)
    ierr = 2;

  // Summarize test results  
  if (ierr == 0)
    utils.out() << "Test passed!" << std::endl;
  else 
    utils.out() << "Test failed!" << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ierr;
}



//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
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
//#include "AztecOO.h"

// User's application specific files 
#include "1DfemInterface.H" 

#include "Teuchos_dyn_cast.hpp"

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
  TransientInterface(int my_NumGlobalElements, Epetra_Comm& my_Comm, 
                     double my_xmin = 0.0, double my_xmax = 1.0) :
    Interface(my_NumGlobalElements, my_Comm, my_xmin, my_xmax),
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
  // Remove shadow warning:
  bool evaluate(NOX::Epetra::Interface::Required::FillType, const Epetra_Vector*, Epetra_Vector*, Epetra_RowMatrix*) {return false;}
  // New function:
  bool evaluate(NOX::Epetra::Interface::Required::FillType flag, 
                const Epetra_Vector* soln, 
                const Epetra_Vector* solndot, // 09/21/05 tscoffe
                double alpha, // 09/21/05 tscoffe
                double beta, // 09/21/05 tscoffe
  		Epetra_Vector* tmp_rhs, 
  		Epetra_RowMatrix* tmp_matrix)
  {

    Epetra_CrsMatrix* jacout = NULL;
    
    //Determine what to fill (F or Jacobian)
    bool fillF = false;
    bool fillMatrix = false;
    if (tmp_rhs != 0) {
      fillF = true;
      rhs = tmp_rhs;
    }
    else {
      fillMatrix = true;
      jacout = &Teuchos::dyn_cast<Epetra_CrsMatrix>(*tmp_matrix);
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
    Epetra_Vector udot(*OverlapMap); // 09/21/05 tscoffe
  
    // Export Solution to Overlap vector
    u.Import(*soln, *Importer, Insert);
    uold.Import(*oldSolution, *Importer, Insert);
    xvec.Import(*xptr, *Importer, Insert);
    udot.Import(*solndot, *Importer, Insert); // 09/21/05 tscoffe
  
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
    double uudot[2]; // 09/21/05 tscoffe
    Basis basis;
  
    // Zero out the objects that will be filled
    if (fillF) 
      rhs->PutScalar(0.0);
    if (fillMatrix) 
      jacout->PutScalar(0.0);
  
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
        uudot[0] = udot[ne]; // 09/21/05 tscoffe
        uudot[1] = udot[ne+1]; // 09/21/05 tscoffe
        // Calculate the basis function at the gauss point
        // basis.computeBasis(gp, xx, uu, uuold);
        basis.computeBasis(gp, xx, uu, uudot); // 09/21/05 tscoffe
  	            
        // Loop over Nodes in Element
        for (int i=0; i< 2; i++) {
  	row=OverlapMap->GID(ne+i);
  	//printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
  	//     MyPID, row, ne+i,StandardMap.MyGID(row));
  	if (StandardMap->MyGID(row)) {
  	  if (fillF) {
            (*rhs)[StandardMap->LID(OverlapMap->GID(ne+i))]+=
              +basis.wt*basis.dx*(
             // (basis.uu-basis.uuold)/dt*basis.phi[i] +   // 09/21/05 tscoffe
                basis.uuold*basis.phi[i] +                 // 09/21/05 tscoffe
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
             // basis.phi[j]/dt*basis.phi[i] +    // 09/21/05 tscoffe
                basis.phi[j]*alpha*basis.phi[i] + // 09/21/05 tscoffe
          beta*(                                  // 09/21/05 tscoffe
	        (1.0/(basis.dx*basis.dx))*
                basis.dphide[j]*basis.dphide[i] -
                8.0/factor/factor*
                (2.0*basis.uu-3.0*basis.uu*basis.uu)*
                basis.phi[j]*basis.phi[i])); // 09/21/05 tscoffe
  	      ierr=jacout->SumIntoGlobalValues(row, 1, &jac, &column);
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
        jacout->ReplaceGlobalValues(0, 1, &jac, &column);
        column=1;
        jac=0.0;
        jacout->ReplaceGlobalValues(0, 1, &jac, &column);
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
      jacout->ReplaceGlobalValues(row, 1, &jac, &column);
      column--;
      jac=0.0;
      jacout->ReplaceGlobalValues(row, 1, &jac, &column);
    }
  }

    // Sync up processors to be safe
    Comm->Barrier();
   
    if (fillMatrix)
      jacout->FillComplete();
  
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


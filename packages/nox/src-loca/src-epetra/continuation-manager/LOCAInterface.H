#ifndef LOCA_INTERFACE_H
#define LOCA_INTERFACE_H

// Trilinos headers
#include <Epetra_Vector.h>
#include <LOCA_Epetra.H>

// ParaCont headers
#include "ProblemLOCAPrototype.H"
#include "ContinuationManager.H"

class  LOCAInterface : 
  public LOCA::Epetra::Interface::Required,
//  public LOCA::Epetra::Interface::TimeDependent,
  public NOX::Epetra::Interface::Jacobian {

    public:

      // Constructor
      LOCAInterface( Teuchos::RCP <ProblemLOCAPrototype> & aProblem ,
	  Teuchos::RCP <ContinuationManager> aContinuationManager);

      // Destructor
      ~LOCAInterface();
      
      // Compute and return F
      bool computeF(const Epetra_Vector& x, Epetra_Vector& f, 
	  NOX::Epetra::Interface::Required::FillType F);
    
      // Compute an explicit Jacobian
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      // Set a parameter in the user's code.
      void setParameters(const LOCA::ParameterVector& params);

      // Print the solution
      void printSolution (const Epetra_Vector &x, const double conParam);

//      // Compute the Shifted matrix needed for linear stability
//      bool computeShiftedMatrix (double alpha, double beta, 
//                                 const Epetra_Vector &x, Epetra_Operator &A);

    private:

      // The continuation manager
      Teuchos::RCP <ContinuationManager> continuationManager;

      // The Problem
      Teuchos::RCP <ProblemLOCAPrototype> problem;

};

#endif

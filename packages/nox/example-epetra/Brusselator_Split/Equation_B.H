// Finite Element Problem Class
/* Provides function (F) and Jacobian evaluations for the following equation 
 * via a 1D linear finite element discretization with Epetra objects.
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */

#ifndef _NOX_EXAMPLE_EPETRA_BRUSSELATOR_EQUATION_B_H
#define _NOX_EXAMPLE_EPETRA_BRUSSELATOR_EQUATION_B_H

#include "GenericEpetraProblem.H" // Base class

// Single PDE Equation Problem Class for Brusselator system

class Equation_B : public GenericEpetraProblem { 

public:

  // Constructor
  Equation_B(int NumGlobalUnknowns, Epetra_Comm& Comm);
  
  // Destructor
  ~Equation_B();

  // Reset problem for next parameter (time) step.
  // For now, this simply updates oldsoln with the given Epetra_Vector
  void reset(const Epetra_Vector& x);

  void reset(); // Empty reset method needed here to prevent hiding base
                // class method.  

  // Set initial condition for solution vector
  void initializeSolution();

  // Evaluates the function (F) and/or the Jacobian using the solution 
  // values in solnVector.
  bool evaluate(FillType f, const Epetra_Vector *solnVector, 
		Epetra_Vector *rhsVector, Epetra_RowMatrix *matrix,
                NOX::Epetra::Interface::FillType fill=NOX::Epetra::Interface::F);

  // Accesor function for time step
  double getdt();

  // Return a reference to the Epetra_Vector with the old solution
  Epetra_Vector& getOldSoln();


protected:

  // inserts the global column indices into the Graph
  void generateGraph();

protected:

  double xmin;
  double xmax;
  double dx;
  double dt;

  FillType flag;
  Epetra_Import *ColumnToOverlapImporter;
  Epetra_Vector *oldSolution;
  Epetra_Vector *rhs;

};
#endif





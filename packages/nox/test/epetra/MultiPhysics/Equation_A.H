// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#ifndef _NOX_EXAMPLE_EPETRA_BRUSSELATOR_EQUATION_A_H
#define _NOX_EXAMPLE_EPETRA_BRUSSELATOR_EQUATION_A_H

#include "GenericEpetraProblem.H" // Base class

// Single PDE Equation Problem Class for Brusselator system

class Equation_A : public GenericEpetraProblem {

public:

  // Constructor
  Equation_A(Epetra_Comm& Comm, int NumGlobalUnknowns = 0,
                       std::string name = "");

  // Destructor
  ~Equation_A();

  // Reset problem for next parameter (time) step.
  // For now, this simply updates oldsoln with the given Epetra_Vector
  void reset(const Epetra_Vector& x);

  void reset(); // Empty reset method needed here to prevent hiding base
                // class method.

  // Set initial condition for solution vector
  void initializeSolution();

  // Initialize based on registrations
  void initialize();

  // Evaluates the function (F) and/or the Jacobian using the solution
  // values in solnVector.
  bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
                const Epetra_Vector *solnVector,
        Epetra_Vector *rhsVector);

  // Accesor function for time step
  double getdt() const { return dt; }

  // Accesor function for setting time step
  void setdt( double dt_ ) { dt = dt_;}

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

  Epetra_Import *ColumnToOverlapImporter;
  Epetra_Vector *oldSolution;
  Epetra_Vector *rhs;

  int id_spec ;
  int id_vel  ;

  bool useConvection;
};
#endif

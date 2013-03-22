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

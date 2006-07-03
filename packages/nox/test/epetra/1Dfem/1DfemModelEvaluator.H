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

#ifndef TEST_EPETRA_1DFEM_MODELEVALUATOR_EPETRA_H
#define TEST_EPETRA_1DFEM_MODELEVALUATOR_EPETRA_H

#include "EpetraExt_ModelEvaluator.h" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#include "Teuchos_RefCountPtr.hpp"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map; 
class Epetra_Vector;
class Epetra_Import;
class Epetra_Operator;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
namespace NOX {
  namespace Parameter {
    class List;
  }
}

// Finite Element Problem Class
class ModelEvaluatorInterface : public EpetraExt::ModelEvaluator { 

public:

  // Constructor
  ModelEvaluatorInterface(int NumGlobalElements, Epetra_Comm& Comm, 
			  double xmin = 0.0, double xmax = 1.0 );
  
  // Destructor
  virtual ~ModelEvaluatorInterface();

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

  // Set PDE factor
  bool setPDEfactor(double value) { factor = value; return true; };

protected:

  // Evaluates the function (F) and/or the Jacobian using the solution 
  // values in solnVector.
  virtual bool evaluate(const Epetra_Vector *solnVector, 
			Epetra_Vector *rhsVector, 
			Epetra_RowMatrix *matrix) const;

  // Inserts the global column indices into the Graph
  virtual bool createGraph();

  // Set initial condition for solution vector
  virtual bool initializeSoln();

protected:

  enum TargetObjectType {FUNCTION, JACOBIAN};

  int NumGlobalElements;  // Total Number of elements
  int NumMyElements;      // Number of elements owned by this process
  int MyPID;              // Process number
  int NumProc;            // Total number of processes

  double xmin;
  double xmax;
  double factor; 	  // PDE factor

  Epetra_Comm* Comm;
  Teuchos::RefCountPtr<Epetra_Map> StandardMap; 
  Teuchos::RefCountPtr<Epetra_Map> OverlapMap;
  Epetra_Import* Importer;
  Teuchos::RefCountPtr<Epetra_Vector> initialSolution;
  mutable Epetra_Vector* rhs;
  Epetra_CrsGraph* Graph;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> jacobian;
  Teuchos::RefCountPtr<Epetra_Vector> xptr;

};

//==================================================================
// Finite Element Basis Object
class Basis {

 public:
  // Constructor
  Basis();

  // Destructor
  ~Basis();

  // Calculates the values of u and x at the specified gauss point
  void computeBasis(int gp, double *x, double *u, double *uold = 0);

 public:
  // Variables that are calculated at the gauss point
  double *phi, *dphide; 
  double uu, xx, duu, eta, wt;
  double dx;
  // These are only needed for transient
  double uuold, duuold;
};

#endif





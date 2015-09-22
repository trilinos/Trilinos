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

#ifndef TEST_EPETRA_1DFEM_MODELEVALUATOR_EPETRA_H
#define TEST_EPETRA_1DFEM_MODELEVALUATOR_EPETRA_H

#include "EpetraExt_ModelEvaluator.h" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#include "Teuchos_RCP.hpp"

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
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
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
  double factor;       // PDE factor

  Epetra_Comm* Comm;
  Teuchos::RCP<Epetra_Map> StandardMap;
  Teuchos::RCP<Epetra_Map> OverlapMap;
  Epetra_Import* Importer;
  Teuchos::RCP<Epetra_Vector> initialSolution;
  mutable Epetra_Vector* rhs;
  Epetra_CrsGraph* Graph;
  Teuchos::RCP<Epetra_CrsMatrix> jacobian;
  Teuchos::RCP<Epetra_Vector> xptr;

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





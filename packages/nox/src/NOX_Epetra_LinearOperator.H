// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_LINEAR_OPERATOR_H
#define _NOX_LINEAR_OPERATOR_H

#include "NOX_Parameter_List.H"
#include "NOX_Epetra_SharedJacobian.H"
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

namespace NOX {
namespace Epetra {


//! Controls the linear solver.
/*! This creates the SharedJacobian and, controls the options for 
 * the linear solver:
 *
 * 1. computeJacobian() calls - If a "Matrix-Free" or "Finite Difference" 
 * jacobian is used, a special computeJacobian() associated with that 
 * matrix must be called instead of the NOX::Epetra::Interface method.
 *
 * 2. All preconditioning options.  NOTE: this has not been fuly 
 * implemented until we get more functionality from 
 * Trilinos (promised by 2/14/0).
 */
class LinearOperator : public Interface {

public:

  //! Constructor with Jacobian Matrix only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_RowMatrix must be supplied for 
   * the Jcobian even if in Matrix-Free mode. */
  LinearOperator(const Parameter::List& params, Interface& i, 
		 Epetra_RowMatrix& J);

  //! Constructor with  Jacobian (J) and Preconditioner (M)
  LinearOperator(const Parameter::List& params, Interface& i, 
		 Epetra_RowMatrix& J, Epetra_RowMatrix& M);

  //! Destructor
  ~LinearOperator();

  //! Solve the linear system in the corresponding Group g. 
  //! Returns true if solved, false if linear solve failed.  
  bool solveLinearSystem(const Parameter::List& params, Group* g);

  //! Set any required Aztec options.  Not all options are present, will be filled in as needed.
  void setAztecOptions(const Parameter::List& params, AztecOO& aztec);

  //! Return the SharedJacobian created by LinearOperator.
  SharedJacobian& getSharedJacobian();

  // Derived from NOX::Epetra::Interface
  bool computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS);
  
  // Derived from NOX::Epetra::Interface
  bool computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac);

  // Derived from NOX::Epetra::Interface
  bool computePreconditioner(Epetra_RowMatrix& M);

  // Derived from NOX::Epetra::Interface
  bool preconditionVector(Epetra_Vector& y);

protected:

  //! Returns the type of Jacobian matrix that is passed into the 
  //! LinearOperator constructors.
  string getJacobianType();

protected:

  //! Determines the Jacobian matrix implementation 
  /*!
    <li> Key: "Jacobian Matrix Type"
    <li> Options:
    <ul> 
    <li> "User Supplied" - user supplied
    <li> "NOX" - created by NOX using Finite Differencing
    <li> "Matrix-Free"  - approximated in NOX using directional derivative
    </ul>
  */ 
  string jacType;
  
  //! Determines the preconditioning matrix implementation 
  /*! 
   <li> Key: "Preconditioning Matrix Type"
   <li> Options:
   <ul>
   <li> "None" - no precondtioning matrix is explicitly defined
   <li> "Use Jacobian" - linear solver library computes preconditioner based on the Jacobian
   <li> "User Supplied Matrix" - linear solver library computes preconditionerbased on supplied matrix
   <li> "User Supplied Action on Vector" - user computes preconditioner and applies it to vector 
   </ul>
  */ 
  string precType;
  
  //! Choice of Aztec preconditioner
  /*! 
    <li> Key: "Preconditioner"
    <li> Options:
    <ul>
    <li> "None"
    <li> "ilut"
    <li> "Polynomial"
    <li> "User" 
    </ul>
  */ 
  string preconditioner;

  //! User supplied NOX::Epetra::Interface
  Interface& userInterface;
  
  //! Jacobian matrix
  Epetra_RowMatrix& jacobian;

  //! Preconditioning matrix
  Epetra_RowMatrix* prec;

  //! Pointer to the SharedJacobian
  SharedJacobian* sharedJacobian;

};
}
}
#endif

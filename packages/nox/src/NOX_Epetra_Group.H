// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_GROUP_H
#define NOX_EPETRA_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Epetra_Vector.H"	// class data element
#include "NOX_Common.H" // class data element (string)

// Forward declares
namespace NOX {
namespace Epetra {
class Interface;
class SharedJacobian;
class LinearOperator;
}
namespace Parameter {
class List;
}
}
class Epetra_Vector;
class Epetra_RowMatrix;
class AztecOO;

namespace NOX {
namespace Epetra {

//! Concrete implementation of NOX::Abstract::Group for %Epetra.
class Group : public Abstract::Group {

public:
  //! Constructor with Jacobian Matrix only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even if in Matrix-Free mode. */
  Group(const Parameter::List& params, Interface& i, 
 	Epetra_Vector& x, Epetra_RowMatrix& J);

  //! Constructor with  Jacobian (J) and Preconditioner (M)
  Group(const Parameter::List& params, Interface& i, 
 	Epetra_Vector& x, Epetra_RowMatrix& J, Epetra_RowMatrix& M);


  /*! \brief Copy constructor. If type is DeepCopy, takes ownership of
    valid shared Jacobian. */
  Group(const Group& source, CopyType type = DeepCopy);

  //! Destructor.
  virtual ~Group();

  virtual Abstract::Group& operator=(const Group& source);
  //! See above.
  virtual Abstract::Group& operator=(const Abstract::Group& source);

  /** @name "Compute" functions. */
  //@{

  virtual bool computeX(const Group& grp, 
			const Vector& d, 
			double step);
  //! See above.
  virtual bool computeX(const Abstract::Group& grp,
			const Abstract::Vector& d, 
			double step);

  virtual bool computeRHS();

  virtual bool computeJacobian();

  virtual bool computeGrad();

  virtual bool computeNewton(Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  virtual bool applyJacobian(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const;

  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  virtual bool applyJacobianDiagonalInverse(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  virtual bool isRHS() const;
  virtual bool isJacobian() const;
  virtual bool isGrad() const;
  virtual bool isNewton() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const Abstract::Vector& getX() const;

  virtual const Abstract::Vector& getRHS() const;
  
  virtual double getNormRHS() const;

  virtual const Abstract::Vector& getGrad() const;

  virtual const Abstract::Vector& getNewton() const;

  //@}

  virtual Abstract::Group* clone(CopyType type = DeepCopy) const;

  //! Return the SharedJacobian created by LinearOperator.
  SharedJacobian& getSharedJacobian();

  //! Return the userInterface owned by LinearOperator.
  Interface& getUserInterface();

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

  //! Set any required Aztec options.  Not all options are present, will be filled in as needed.
  void setAztecOptions(const Parameter::List& params, AztecOO& aztec);

  //! Returns the type of Jacobian matrix that is passed into the 
  //! LinearOperator constructors.
  string getJacobianType();

protected:

  /** @name Vectors */
  //@{
  //! Solution vector.
  Vector xVector;
  //! Right-hand-side vector (function evaluation).
  Vector RHSVector;
  //! Gradient vector (steepest descent vector).
  Vector gradVector;
  //! Newton direction vector.
  Vector NewtonVector;
  //@}

  //! An extra temporary vector, only allocated if needed.
  mutable Epetra_Vector* tmpVectorPtr;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidRHS;
  bool isValidJacobian;
  bool isValidGrad;
  bool isValidNewton;
  //@}
  
  //! Norm of RHS
  double normRHS;

  //! Pointer to shared Jacobian matrix 
  SharedJacobian* sharedJacobianPtr;

  //! Reference to shared Jacobian matrix 
  SharedJacobian& sharedJacobian;

  //! Reference to the user supplied interface functions
  Interface& userInterface;

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

};

} // namespace Epetra
} // namespace NOX


#endif

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
}
namespace Parameter {
class List;
}
}
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;
class AztecOO;

namespace NOX {
namespace Epetra {

//! Concrete implementation of NOX::Abstract::Group for %Trilinos/%Epetra.
/*! This group is set up to use the AztecOO Linear %Solver.
  The user has multiple options for forming the Jacobian operator and preconditioning.  
  
  A linear solve using this Group class requires that the user supply at a minimum an object for the Jacobian operator derived from the pure virtual Epetra_Operator class.  This does not have to be an explicit matrix since Newton-Krylov methods only require the action of the Jacobian on a vector (Jy).  The user can pass in their own Epetra_Operator or use one of the predefined Epetra_Operators provided by %Trilinos/%Epetra (Epetra_CreMatrix Epetra_VbrMatrix) or NOX (NOX::Epetra::MatrixFree or NOX::Epetra::FiniteDifference).  Certain functions in the NOX::Epetra::Interface object are required to be supplied by the user depending on the type of Jacobian used: 
  <ul>
  <li> User Created Epetra_Operator or %Trilinos/%Epetra - The user must supply the operator to the group constructor.  This could  be any Epetra_RowMatrix or a user created class.  To be able to use this object, NOX requires the user to implement the following functions in the NOX::Epetra::Interface object: 
  <ul>
  <li> computeRHS()
  <li> computeJacobian() 
  </ul>

  <li> NOX::Epetra::MatrixFree - This object represents a Matrix-Free implementation an thus requires minimal memory since no Jacobian is stored.  The user need only supply the computeRHS() function for the NOX::Epetra::Interface class.

  <li> NOX::Epetra::FiniteDifference - This object uses finite differencing of the RHS to compute and store the entire Jacobian as an Epetra_CrsMatirx.  This is not optimied and can be very time consuming since an Interface::computeRHS() call is made for each unknown. This was meant primarily as a validation tool to ensure Jacobian fills are correct. The user need only supply the computeRHS() function for the NOX::Epetra::Interface class. 
  </ul>
 
Preconditioning also has multiple options that may depend on the Jacobian. A flag in the Linear Solver sublist called "Preconditioning" must be set to turn on any preconditioning of the linear solver.  A preconditioner object is an Epetra_Operator derived class just like th Jacobian.  If the AztecOO preconditioners are to be used, they require some knowledge of the matrix and thus at a minimum the preconditioner object must also be an Epetra_RowMatrix. The valid options are listed below:
  <ul>
  <li> "None" - no preconditioning.  This is the default setting.
  <li> "AztecOO: Jacobian Matrix" - An AztecOO preconditioner will be used.  The preconditioning will be based on the Jacobian.  This requires that the Jacobian operator MUST BE an Epetra_RowMatrix (NOX will test for this and throw an error if it is not).  This means a NOX::Epetra::MatrixFree Jacobian can NOT be used with this option for preconditioning.
  <li> "AztecOO: User RowMatrix" - An AztecOO preconditioner will be used.  The preconditioning will be based on an Epetra_Operator that is separate from the Jacobian.  The constructor with separate Jacobian and preconditioner Epetra_Operators must be used. Since AztecOO preconditioners require matrix information, the preconditioning object (M) must be an Epetra_RowMatix (NOX will test for this and throw an error if it is not). This option requires the user to implement the function computePrecMatrix() in the NOX::Epetra::Interface derived class.
  <li> "User Supplied Preconditioner" - The user supplies an Epetra_Operator that will compute and apply the preconditioner to a vector that is passed in. This option requires the user to implement the function computePreconditioner() in the NOX::Epetra::Interface derived class.   
  </ul>

If an AztecOO preconditioner is used, the type of preconditioner can be set with the flag "Aztec Preconditioner" in the "Linear Solver" sublist (defaults to ilu). Currently supported preconditioners and their corresponding parameters that can be set in the "Linear Solver" sublist are shown below (See the Aztec 2.0 manual for more information):
  <ul>
  <li> "ilu" 
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Graph Fill" - defaults to 0.
  </ul> 
  <li> "ilut"
  <ul>
  <li> "Overlap" - defaults to 0.
  <li> "Fill Factor" - defaults to 1.0.
  <li> "Drop Tolerance" - defaults to 1.0e-12
  </ul>
  <li> "Polynomial" - Neumann polynomial
  <ul>
  <li> "Polynomial Order" - defaults to 3.
  </ul>
  </ul>

Addtional Parameters that are set in the "Linear Solver" sublist:
  <ul> 
  <li> "Tolerance" - tolerance used by AztecOO to determine if a linear solve has converged.
  <li> "Max Iterations" - maximum number of linear iterations in the 
  <li> "Output Frequency" - number of linear solve iterations between output of the linear solve residual. Takes an integer, AZ_none, or AZ_all as a value.
  </ul>


 */

class Group : public Abstract::Group {

public:
  //! Constructor with Jacobian Operator only.  
  /*! Either there is no preconditioning or the Jacobian will be used
   * for preconditioning.  An Epetra_Operator must be supplied for 
   * the Jacobian even if in Matrix-Free mode. linearSolverParams is 
   * the "Linear Solver" sublist of parameter list.
   */
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J);

  //! Constructor with  Jacobian (J) and a separate Preconditioner (M).  
  //! linearSolverParams is the "Linear Solver" sublist of parameter list.
  Group(const Parameter::List& linearSolverParams, Interface& i, 
 	Epetra_Vector& x, Epetra_Operator& J, Epetra_Operator& M);


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

  //! Return the SharedJacobian.
  virtual SharedJacobian& getSharedJacobian();

  //! Return the userInterface.
  virtual Interface& getUserInterface();

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

  /** @name "Linear Solver" functions 
   */
  //@{
  //! Set any required Aztec options.  NOTE: This does not handle all aztec options.  They will be added as needed. 
  virtual void setAztecOptions(const Parameter::List& params, AztecOO& aztec);
  
  //! Sets the linear solver options and checks to make sure the supplied objects are valid for those options.  
  virtual bool setLinearSolver(const Parameter::List& params);
  
  //! Computes the preconditioner based on the current solution vector and sets the correct preconditioning option (Operator or RowMatrix) in the AztecOO object.
  virtual bool computePreconditioner(AztecOO& aztec);
  
  //! Returns the type of Jacobian matrix that is passed into the 
  //! group constructors.
  virtual string getJacobianType();

  //! Returns the type of Preconditioning matrix that is passed into the 
  //! group constructors.
  virtual string getPrecType();
  //@}

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

  //! Flag that tells the code how the Jacobian operator is implemented.
  /*! This flag is set in the constructor by calling the function getJacobianType().
    <li> Options:
    <ul> 
    <li> "User Supplied" - user supplied Epetra_Operator or Epetra_RowMatrix.
    <li> "Finite Difference" - the Jacobian is a NOX::Epetra::FiniteDifference object (that is really an Epetra_CrsMatrix created using using Finite Differencing).
    <li> "Matrix-Free" - No explicit Jacobian matrix is created.  The Jacobian is a NOX::Epetra::MatrixFree object.  The matrix-vector products required for Newton-Krylov iterative solves are approximated in NOX using directional derivatives.
    </ul>
  */ 
  string jacType;
  
  //! Flag that tells the code how the Preconditioning operator is implemented.
  /*! This flag is set in the constructor by calling the function getPrecType().
    <li> Options:
    <ul> 
    <li> "User Supplied" - user supplied Epetra_Operator or Epetra_RowMatrix.
    <li> "Finite Difference" - NOX::Epetra::FiniteDifference object (that is really an Epetra_CrsMatrix created using using Finite Differencing).
    </ul>
  */ 
  string precType;
  
  //! Determines how the preconditioning is handled. 
  /*! This variable is set in the constructors from the "Preconditioning" parameter in the "Linear Solver" sublist.
   <li> Key in "Linear Solver" sublist: "Preconditioning"
   <li> Valid Options:
   <ul>
   <li> "None" - no precondtioning is used.
   <li> "AztecOO: Jacobian Matrix" - AztecOO preconditioner is used.  This requires an explicit Jacobian Matrix.  Therefore the Jacobian matrix must be an Epetra_RowMatrix or a NOX_FiniteDifference object.
   <li> "AztecOO: User RowMatrix" - AztecOO preconditioner is used. This requires an explicit Epetra_RowMatrix that is NOT the same as the Jacobian.  This was specifically written for "Matrix-Free" mode when there is no explicit Matrix for the Jacobian.  
   <li> "User Supplied Preconditioner" - user supplied routine computes preconditioner and applies it to vector.  The user must implement the routine as an Epetra_Operator class such that the Apply() function preconditions the incoming vector. 
   </ul>
  */ 
  string precOption;
  
  //! Choice of Aztec preconditioner.  
  /*! This variable is set in the constructors from the "Aztec Preconditioner" parameter in the "Linear Solver" sublist.
    <li> Key in "Linear Solver" sublist: "Aztec Preconditioner"
    <li> Valid Options:
    <ul>
    <li> "ilu"
    <li> "ilut"
    <li> "Polynomial" 
    </ul>
  */ 
  string aztecPrec;

};

} // namespace Epetra
} // namespace NOX


#endif

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

#ifndef NOX_PETSC_GROUP_H
#define NOX_PETSC_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Petsc_Vector.H"	// class data element

#include "petscsnes.h"	
#include "petscsles.h"	

// Forward declares
namespace NOX {
namespace Petsc {
class Interface;
class SharedJacobian;
}
namespace Parameter {
class List;
}
}

namespace NOX {
namespace Petsc {

//! Concrete implementation of NOX::Abstract::Group for %Petsc.
class Group : public Abstract::Group {

public:

  //! Constructor.
  Group(const Parameter::List& params, Interface& i,
             Vec& x, Mat& J, string jType);

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

  virtual bool setX(const Vector& y);

  //! See above
  virtual bool setX(const Abstract::Vector& y);

  virtual bool computeX(const Group& grp, 
			const Vector& d, 
			double step);
  //! See above.
  virtual bool computeX(const Abstract::Group& grp,
			const Abstract::Vector& d, 
			double step);

  virtual bool computeF();

  virtual bool computeJacobian();

  virtual bool computeGradient();

  virtual bool computeNewton(Parameter::List& params);

  virtual bool computePrecMatrix();

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
  
  /*! 
      \brief Preconditions a vector using one of the available methods in %Petsc
 
      Currently, the preconditioning method corresponds to whatever
      is available in the user's %Petsc installation and is specified
      in the local \b .petscrc file.  For example, a direct solve can
      be accomplished by specifying \b -pc_type \b lu (serial only),
      whereas diagonal scaling can be achieved via \b -pc_type \b jacobi.
  */ 
  virtual bool applyPrecMatrixInverse(const Vector& input, Vector& result) const;
  //! See above
  virtual bool applyPrecMatrixInverse(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  virtual bool isF() const;
  virtual bool isJacobian() const;
  virtual bool isGradient() const;
  virtual bool isNewton() const;
  virtual bool isPrecMatrix() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  virtual const Abstract::Vector& getX() const;

  virtual const Abstract::Vector& getF() const;
  
  virtual double getNormF() const;

  virtual const Abstract::Vector& getGradient() const;

  virtual const Abstract::Vector& getNewton() const;

  //@}

  virtual Abstract::Group* clone(CopyType type = DeepCopy) const;

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

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
  mutable Vec* tmpVectorPtr;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidRHS;
  bool isValidJacobian;
  bool isValidGrad;
  bool isValidNewton;
  bool isValidPrecMatrix;
  //@}
  
  //! Norm of RHS
  double normRHS;

  //! Pointer to shared Jacobian matrix
  SharedJacobian* sharedJacobianPtr;

  //! Reference to shared Jacobian matrix 
  SharedJacobian& sharedJacobian;

  //! Specification of Jacobian Type
  string jacType;

  //! Reference to the user supplied interface functions
  Interface& userInterface;

};

} // namespace Petsc
} // namespace NOX



#endif

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

#ifndef NOX_EXAMPLE_GROUP_H
#define NOX_EXAMPLE_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "NOX_Common.H" // class data element (string)
#include "NOX_Example_Vector.H"	// class data element
#include "NOX_Example_Matrix.H"
#include "NOX_Example_Interface.H"

// Forward declares
namespace NOX {
namespace Parameter {
class List;
}
}

namespace NOX {
namespace Example {

//! A simple example of a group structure, based on LAPACK.
class Group : public Abstract::Group {

public:

  //! Constructor
  Group(Interface& i);

  //! Copy constructor
  Group(const Group& source, CopyType type = DeepCopy);

  //! Destructor.
  ~Group();

  Abstract::Group& operator=(const Group& source);

  //! See above.
  Abstract::Group& operator=(const Abstract::Group& source);

  /** @name "Compute" functions. */
  //@{

  bool setX(const Vector& y);

  //! See above
  bool setX(const Abstract::Vector& y);

  bool computeX(const Group& grp, const Vector& d, double step);
  //! See above.
  bool computeX(const Abstract::Group& grp, const Abstract::Vector& d, double step);

  bool computeF();

  bool computeJacobian();

  bool computeGradient();

  bool computeNewton(Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  bool applyJacobian(const Vector& input, Vector& result) const;

  //! See above
  bool applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const;

  bool applyJacobianTranspose(const Vector& input, Vector& result) const;

  //! See above
  bool applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  bool applyJacobianDiagonalInverse(const Vector& input, Vector& result) const;

  //! See above
  bool applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const;
  
  bool preconditionVector(const Vector& input, Vector& result) const;

  //! See above
  bool preconditionVector(const Abstract::Vector& input, Abstract::Vector& result) const;

  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  bool isF() const;
  bool isJacobian() const;
  bool isGradient() const;
  bool isNewton() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  const Abstract::Vector& getX() const;

  const Abstract::Vector& getF() const;
  
  double getNormF() const;

  const Abstract::Vector& getGradient() const;

  const Abstract::Vector& getNewton() const;

  //@}

  Abstract::Group* clone(CopyType type = DeepCopy) const;

  //! Print out the group
  void print() const;

private:

  //! resets the isValid flags to false
  void resetIsValid();

private:

  /** @name Vectors */
  //@{
  //! Solution vector.
  Vector xVector;
  //! Right-hand-side vector (function evaluation).
  Vector fVector;
  //! Newton direction vector.
  Vector newtonVector;
  //! Gradient vector (steepest descent vector).
  Vector gradientVector;
  //@}


  //! Jacobian Matrix
  Matrix jacobianMatrix;


  //! Problem interface
  Interface& problemInterface;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidF;
  bool isValidJacobian;
  bool isValidGradient;
  bool isValidNewton;
  //@}
  
  //! Norm of F
  double normF;


};

} // namespace Example
} // namespace NOX


#endif

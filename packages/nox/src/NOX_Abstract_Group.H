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

#ifndef NOX_ABSTRACT_GROUP_H
#define NOX_ABSTRACT_GROUP_H

#include "NOX_Abstract_Vector.H" // for NOX::CopyType
#include "NOX_Common.H" // for string 

namespace NOX {
namespace Parameter {
class List;
}
}

namespace NOX { 
namespace Abstract { 

/*! 
  \brief %NOX pure abstract interface to a "group"; i.e., a solution
  vector and the corresponding RHS, Jacobian, gradient, and Newton
  vectors, as well as any shared objects.
 
  This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.
*/

class Group {

public:
  
  //! Constructor.
  Group() {};

  //! Destructor.
  virtual ~Group() {};
  
  /*! 
    \brief Copies the values of all vectors and any other data in
    source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual Group& operator=(const Group& source) = 0;
      
  //@{ \name "Compute" functions.

  //! Compute and return solution vector, x, where this.x = grp.x() + step * d.
  virtual bool computeX(const Group& grp, const Vector& d, double step) = 0;

  //! Compute and return RHS. (Usually also computes and stores norm of RHS.)
  virtual bool computeRHS() = 0;

  //! Compute Jacobian.
  virtual bool computeJacobian() = 0;

  //! Compute and return gradient.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual bool computeGrad() = 0;

  //! Compute and return Newton direction, using parameters for nonlinear solve.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual bool computeNewton(Parameter::List& params) = 0;

  //@}

  /** @name Jacobian operations.
   
    Operations using the Jacobian matrix. These may not be defined in
    matrix-free scenarios.
  */

  //@{
  
  /*! 
    \brief If supported, returns true and calculates 
    result = Jacobian * input.  Otherwise, returns false.  
    Returns false if any errors occur, such as the Jacobian not being
    computed.
  */
  virtual bool applyJacobian(const Vector& input, Vector& result) const = 0;

  /*! 
    \brief If supported, returns true and calculates result =
    Jacobian^T * input.  Otherwise, returns false.  Throws an error if
    the Jacobian has not been computed. 
  */
  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const = 0;
  
  /*!
    \brief Applies the Jacobian Diagonal to the given input vector.
  */
  virtual bool applyJacobianDiagonalInverse(const Vector& input, Vector& result) const = 0;
    

  //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns true
    if the corresponding "compute" function has been called since the
    last update to the solution vector (via instantiation or
    computeX).
  */

  //@{
  
  //! Return true if the RHS is valid.
  virtual bool isRHS() const = 0;
  //! Return true if the Jacobian is valid.
  virtual bool isJacobian() const = 0;
  //! Return true if the gradient is valid.
  virtual bool isGrad() const = 0;
  //! Return true if the Newton direction is valid.
  virtual bool isNewton() const = 0;

  //@}

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const Vector& getX() const = 0;

  //! Return right-hand-side (RHS). 
  virtual const Vector& getRHS() const = 0;

  //! Return 2-norm of RHS.
  virtual double getNormRHS() const = 0;

  //! Return gradient.
  virtual const Vector& getGrad() const = 0;

  //! Return Newton direction.
  virtual const Vector& getNewton() const = 0;

  //@}


  //@{ \name Creating new Groups.

  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is "DeepCopy", then we need to create an exact replica of
    "this". Otherwise, if type is "CopyShape", we need only replicate
    the shape of "this". Returns NULL if clone is not supported.
  */
  virtual Group* clone(CopyType type = DeepCopy) const = 0;

  //@}

};
} // namespace Abstract
} // namespace NOX

#endif

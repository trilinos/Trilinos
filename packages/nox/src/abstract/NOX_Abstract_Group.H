// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_ABSTRACT_GROUP_H
#define _NOX_ABSTRACT_GROUP_H

#include "NOX_Parameter_List.H"	
#include "NOX_Abstract_Vector.H"

namespace NOX { 
namespace Abstract { 

//! 

/*! \brief NOX pure abstract interface to a "group"; i.e., a solution
 *  vector and the corresponding RHS, Jacobian, gradient, and Newton
 *  vectors, as well as any shared objects.
 *
 * This class groups together the objects relating to a particular
 * solution vector. Each time the solution vector is updated via the
 * "computeX" objects, all the previously computed values become
 * invalid. */

class Group {

public:
  
  //! Constructor.
  Group() {};

  //! Destructor.
  virtual ~Group() {};
  
  /*! \brief Copies the values of all vectors and any other data in
   *  source group to this group.  (May invalidate shared data for
   *  source group.) */
  virtual Group& operator=(const Group& source) = 0;
      
  //@{ \name "Compute" functions.

  //! Compute and return solution vector, x, where this.x = grp.x() + step * d.
  virtual const Vector& computeX(const Group& grp, const Vector& d, double step) = 0;

  //! Compute Jacobian.
  virtual void computeJacobian() = 0;

  //! Compute and return RHS. (Usually also computes and stores norm of RHS.)
  virtual const Vector& computeRHS() = 0;

  //! Compute and return gradient.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual const Vector& computeGrad() = 0;

  //! Compute and return Newton direction, using specified parameters for nonlinear solve.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual const Vector& computeNewton(Parameter::List& params) = 0;

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  /*! \brief If supported, returns true and calculates result =
   *  Jacobian * input.  Otherwise, returns false.  Throws an error if
   *  the Jacobian has not been computed. */
  virtual bool applyJacobian(const Vector& input, Vector& result) const = 0;

  /*! \brief If supported, returns true and calculates result =
   *  Jacobian^T * input.  Otherwise, returns false.  Throws an error
   *  if the Jacobian has not been computed. */
  virtual bool applyJacobianTranspose(const Vector& input, Vector& result) const = 0;
  
  //@}

  /** @name "Is" functions.
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */

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
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
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

  /*! \brief Create a new Group of the same derived type as this one
   * by cloning this one, and return a pointer to the new group.  If
   * type is "DeepCopy", then we need to create an exact replica of
   * "this". Otherwise, if type is "CopyShape", we need only replicate
   * the shape of "this". Returns NULL if clone is not supported. */
  virtual Group* clone(CopyType type = DeepCopy) const = 0;

  //@}

};
}
}

#endif

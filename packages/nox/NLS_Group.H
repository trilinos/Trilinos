// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_GROUP_H
#define _NLS_GROUP_H

#include "NLS_ParameterList.H"
#include "NLS_Vector.H"

//! NLSPACK pure abstract interface to a "group"; i.e., a solution
//! vector and the corresponding RHS, Jacobian, gradient, and Newton
//! direction.

/*! This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.  An interesting note is that this group structure does not
  preclude having just one Jacobian matrix - it just means we need to
  keep track of which solution group that Jacobian belongs too. This
  is done in the concrete instantiations of this class. */

class NLS_Group {

protected:
  //! NLS_Group constructor
  NLS_Group() {};

public:

  //! NLS_Group deconstructor
  virtual ~NLS_Group() {};

  //! Copies the values of all vectors in source group to this group
  virtual NLS_Group& copy(const NLS_Group& copyFrom) = 0;
  
  //@{ \name "Compute" functions

  //! Compute and return solution vector, x, where this.x = grp.x() + step * d
  virtual const NLS_Vector& computeX(const NLS_Group& grp, const NLS_Vector& d, double step) = 0;

  //! Compute and return RHS
  virtual const NLS_Vector& computeRHS() = 0;

  //! Returns true if this group is allowed to store a Jacobian matrix.
  /*! Not every group is allowed to own a Jacobian because of its
    size. Typically, only one group at a time will need to own a
    Jacobian, and that will be in a Jacobian-enabled group. Obviously
    a group cannot compute the Gradient or Newton direction unless it
    is Jacobian-enabled. */
  virtual bool isJacobianEnabled() const = 0;

  //! Compute Jacobian (return true if successful)
  virtual void computeJacobian() = 0;

  //! Compute and return gradient 
  //! Throws an error if RHS and Jacobian have not been computed 
  virtual const NLS_Vector& computeGrad() = 0;

  //! Compute and return Newton direction 
  //! Throws an error if RHS and Jacobian have not been computed 
  virtual const NLS_Vector& computeNewton() = 0;

  //! Compute and return Newton direction, using specified parameters for nonlinear solve
  //! Throws an error if RHS and Jacobian have not been computed 
  virtual const NLS_Vector& computeNewton(NLS_ParameterList& params) = 0;

  //! Compute and return Linear RHS vector d.
  virtual double computeLinearRHSNorm() = 0;

  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  virtual bool isRHS() const = 0;
  virtual bool isJacobian() const = 0;
  virtual bool isGrad() const = 0;
  virtual bool isNewton() const = 0;

  //@}

  //@{ \name "Get" functions 

  //! Return solution vector.  This should be protected!  
  //! Getting a pointer to the solution vector could allow users to change
  //! its value without using computeX().  This breaks the isRHS(), 
  //! isJacobian(), isGrad(), and isNewton() functions. 
  virtual const NLS_Vector& getX() const = 0;

  //! Return rhs (throws an error if RHS has not been computed)
  virtual const NLS_Vector& getRHS() const = 0;

  //! Return gradient (throws an error if gradient has not been computed)
  virtual const NLS_Vector& getGrad() const = 0;

  //! Return Newton direction (throws an error if newton direction has not been computed)
  virtual const NLS_Vector& getNewton() const = 0;

  //@}


  //@{ \name Creating new NLS_Groups

  //! Create a \b new NLS_Group of the same derived type as this one, but \b without an enabled Jacobian.
  virtual NLS_Group* newCopy(bool isJacobianEnabled = false) const = 0;

  //@}
};


#endif

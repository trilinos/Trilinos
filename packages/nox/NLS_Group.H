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

  //! Copies the values of all vectors and Jacobian (if it exists) in source group to this group.
  //! (May invalidate Jacobian for source group.)
  virtual NLS_Group& operator=(const NLS_Group& copyFrom) = 0;

  //@{ \name "Compute" functions

  //! Compute and return solution vector, x, where this.x = grp.x() + step * d.
  virtual const NLS_Vector& computeX(const NLS_Group& grp, const NLS_Vector& d, double step) = 0;

  //! Compute and return RHS. (Usually also computes and stores norm of RHS.)
  virtual const NLS_Vector& computeRHS() = 0;

  //! Compute Jacobian.
  virtual void computeJacobian() = 0;

  //! Compute and return gradient.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual const NLS_Vector& computeGrad() = 0;

  //! Compute and return Newton direction, using specified parameters for nonlinear solve.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual const NLS_Vector& computeNewton(NLS_ParameterList& params) = 0;

  //! Compute and return norm of predicted RHS for given update vector.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual double computeNormPredictedRHS(const NLS_Vector& update) = 0;

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

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{ 

  //! Return solution vector.  
  virtual const NLS_Vector& getX() const = 0;

  //! Return right-hand-side (RHS). 
  virtual const NLS_Vector& getRHS() const = 0;

  //! Return 2-norm of RHS.
  virtual double getNormRHS() const = 0;

  //! Return gradient.
  virtual const NLS_Vector& getGrad() const = 0;

  //! Return Newton direction.
  virtual const NLS_Vector& getNewton() const = 0;

  //@}


  //@{ \name Creating new NLS_Groups

  //! Create a \b new NLS_Group of the same derived type as this one.
  //! (May invalidate Jacobian for this group.)
  virtual NLS_Group* newCopy() const = 0;

  //@}
};


#endif

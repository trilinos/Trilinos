// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_CONCRETE_GROUP_H
#define _NLS_CONCRETE_GROUP_H

#include "NLS_Parameter.H"
#include "NLS_Vector.H"
#include "NLS_Group.H"

//! Concrete implementation of the pure virtual NLS_Group.

/*! This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.  An interesting note is that this group structure does not
  preclude having just one Jacobian matrix - it just means we need to
  keep track of which solution group that Jacobian belongs too. This
  is done in the concrete instantiations of this class. */

class NLS_GenericGroup : public NLS_Group {

public:

  //! NLS_Group constructor
  NLS_GenericGroup();

  //! Copy constructor
  NLS_GenericGroup(NLS_Group& copyFrom);

  //! Create a new group where the new solution vector is grp.x() + step * d
  NLS_GenericGroup(NLS_Group& grp, NLS_Vector& d, double step);
 
  //! NLS_Group deconstructor
  virtual ~NLS_GenericGroup();

  //! Copy
  virtual NLS_GenericGroup& operator=(NLS_GenericGroup& copyFrom);

  /** @name "Compute" functions */
  //@{

  //! Compute and return solution vector
  virtual NLS_Vector& computeX(NLS_Group& x, NLS_Vector& d, double step);

  //! Compute and return RHS
  virtual NLS_Vector& computeRHS();

  //! Compute RHS
  virtual void computeJacobian();

  //! Compute and return gradient 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeGrad();

  //! Compute and return Newton direction 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeNewton();

  //! Compute and return Newton direction, using desired accuracy for nonlinear solve
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeNewton(string& name, NLS_Parameter& parameter);

  //@}

  /** @name Checks to see if various objects have been computed. 
   *
   * Returns true if the corresponding "compute" function has been
   * called since the last update to the solution vector (via
   * instantiation or computeX). */
  //@{

  virtual bool isRHS();
  virtual bool isJacobian();
  virtual bool isGrad();
  virtual bool isNewton();

  //@}

  /** @name "Get" function */
  //@{

  //! Return solution vector
  virtual NLS_Vector& getX();

  //! Return rhs (throws an error if RHS has not been computed)
  virtual NLS_Vector& getRHS();

  //! Return gradient (throws an error if gradient has not been computed)
  virtual NLS_Vector& getGrad();

  //! Return Newton direction (throws an error if newton direction has not been computed)
  virtual NLS_Vector& getNewton();

  //@}

  

};


#endif

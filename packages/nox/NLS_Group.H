#ifndef _NLS_GROUP_H
#define _NLS_GROUP_H

#include "NLS_Parameter.H"
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
  NLS_Group();
  
  /** @name Recommended constructors formats */
  //@{
  //! Copy constructor
  NLS_Group(NLS_Group& copyFrom);
  //! Create a new group where the new solution vector is grp.x() + step * d
  NLS_Group(NLS_Group& grp, NLS_Vector& d, double step);
  //@}

public:
  
  //! NLS_Group deconstructor
  virtual ~NLS_Group() = 0;

  //! Copy
  virtual NLS_Group& operator=(NLS_Group& copyFrom) = 0;

  /** @name "Compute" functions */
  //@{

  //! Compute and return solution vector
  virtual NLS_Vector& computeX(NLS_Group& x, NLS_Vector& d, double step) = 0;

  //! Compute and return RHS
  virtual NLS_Vector& computeRHS() = 0;

  //! Compute RHS
  virtual void computeJacobian() = 0;

  //! Compute and return gradient 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeGrad() = 0;

  //! Compute and return Newton direction 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeNewton() = 0;

  //! Compute and return Newton direction, using desired accuracy for nonlinear solve
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual NLS_Vector& computeNewton(string& name, NLS_Parameter& parameter) = 0;

  //@}

  /** @name Checks to see if various objects have been computed. 
   *
   * Returns true if the corresponding "compute" function has been
   * called since the last update to the solution vector (via
   * instantiation or computeX). */
  //@{

  virtual bool isRHS() = 0;
  virtual bool isJacobian() = 0;
  virtual bool isGrad() = 0;
  virtual bool isNewton() = 0;

  //@}

  /** @name "Get" function */
  //@{

  //! Return solution vector.  This should be protected!  
  //! Getting a pointer to the solution vector could allow users to change
  //! its value without using computeX().  This breaks the isRHS(), 
  //! isJacobian(), isGrad(), and isNewton() functions. 
  virtual NLS_Vector& getX() = 0;

  //! Return rhs (throws an error if RHS has not been computed)
  virtual NLS_Vector& getRHS() = 0;

  //! Return gradient (throws an error if gradient has not been computed)
  virtual NLS_Vector& getGrad() = 0;

  //! Return Newton direction (throws an error if newton direction has not been computed)
  virtual NLS_Vector& getNewton() = 0;

  //@}

  

};


#endif

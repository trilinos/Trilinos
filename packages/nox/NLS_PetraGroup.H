// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PETRA_GROUP_H
#define _NLS_PETRA_GROUP_H

#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

//#include "AztecOO.h"
#include "NLS_Parameter.H"
#include "NLS_ParameterList.H"
#include "NLS_Vector.H"
#include "NLS_PetraVector.H"
#include "NLS_Group.H"
#include "NLS_PetraGroupInterface.H"

//! Concrete implementation of the pure virtual NLS_Group.

/*! This class groups together the objects relating to a particular
  solution vector. Each time the solution vector is updated via the
  "computeX" objects, all the previously computed values become
  invalid.  An interesting note is that this group structure does not
  preclude having just one Jacobian matrix - it just means we need to
  keep track of which solution group that Jacobian belongs too. This
  is done in the concrete instantiations of this class. */

class NLS_PetraGroup : public NLS_Group {

public:

  //! Constructor for Jacobian-Enabled group.
  /*! Makes internal xVector a deep copy of x and the internal Jacobian
    a pointer to J. Does not delete J. */
  NLS_PetraGroup(Epetra_Vector& x, Epetra_RowMatrix& J, NLS_PetraGroupInterface& I);

  //! Constructor for Jacobian-Diabled group.
  /*! Makes internal xVector a deep copy of x. A group constructed in
    this manner cannot compute its Jacobian. */
  NLS_PetraGroup(const NLS_PetraGroup& copyFrom);

  //! Deconstructor
  virtual ~NLS_PetraGroup();

  virtual NLS_Group& copy(const NLS_Group& copyFrom);

  //! Copy everything except the Jacobian matrix, if any.
  virtual NLS_Group& copy(const NLS_PetraGroup& copyFrom);

  /** @name "Compute" functions */
  //@{

  //! Compute and return solution vector
  virtual const NLS_Vector& computeX(const NLS_Group& grp, const NLS_Vector& d, double step);

  //! Compute and return solution vector
  virtual const NLS_Vector& computeX(const NLS_PetraGroup& grp, const NLS_PetraVector& d, double step);

  //! Compute and return RHS
  virtual const NLS_Vector& computeRHS();

  //! Compute Jacobian
  virtual void computeJacobian();

  //! Compute and return gradient 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual const NLS_Vector& computeGrad();

  //! Compute and return Newton direction 
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual const NLS_Vector& computeNewton();

  //! Compute and return Newton direction, using desired accuracy for nonlinear solve
  /*! Throws an error if RHS and Jacobian have not been computed */
  virtual const NLS_Vector& computeNewton(NLS_ParameterList& params);

  //! Compute and place linear residual in vector d.
  /*! Compute the RHS based on the local linear model.  The Jacobian and
   *  Newton direction must be alredy evaluated for the solution. This
   *  is typically used for adaptively changing the linear solver tolerance 
   *  for iterative linear solution methods. */  
  virtual double computeLinearRHSNorm();  

  //@}

  //! Returns true if this group can compute a Jacobian
  bool isJacobianEnabled() const;

  /** @name Checks to see if various objects have been computed. 
   *
   * Returns true if the corresponding "compute" function has been
   * called since the last update to the solution vector (via
   * instantiation or computeX). */
  //@{

  virtual bool isRHS() const;
  virtual bool isJacobian() const;
  virtual bool isGrad() const;
  virtual bool isNewton() const;

  //@}

  /** @name "Get" function */
  //@{

  //! Return solution vector
  virtual const NLS_Vector& getX() const;

  //! Return rhs (throws an error if RHS has not been computed)
  virtual const NLS_Vector& getRHS() const;

  //! Return gradient (throws an error if gradient has not been computed)
  virtual const NLS_Vector& getGrad() const;

  //! Return Newton direction (throws an error if newton direction has not been computed)
  virtual const NLS_Vector& getNewton() const;

  //@}

  virtual NLS_Group* newCopy(bool isJacobianEnabled = false) const;

protected:

  //! Disallow Copy operator for now.
  virtual NLS_PetraGroup& operator=(NLS_PetraGroup& copyFrom) {};

  //! resets the isValid flags to false
  void reset();

/** @name Vectors */
//@{
//! Solution vector.
  NLS_PetraVector xVector;
//! Right-hand-side vector (function evaluation).
  NLS_PetraVector RHSVector;
//! Gradient vector (steepest descent vector).
  NLS_PetraVector gradVector;
//! Newton direction vector.
  NLS_PetraVector NewtonVector;
//@}

  //! Pointer to Epetra Jacobian Matrix 
  Epetra_RowMatrix* Jac;

  //! Pointer to user's code for filling the RHS (function) and Jacobian
  NLS_PetraGroupInterface& Interface;

/** @name IsValid flags 
 *  Tells if the vector is out of date with respect to the current 
 *  solution vector.
 *  false = vector is now out of date.
 *  true = vector is up to date with current solution. */
//@{
bool isValidRHS;
bool isValidJacobian;
bool isValidGrad;
bool isValidNewton;
//@}

};


#endif

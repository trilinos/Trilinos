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

#include "NLS_Group.H"		// base class
#include "NLS_ParameterList.H"
#include "NLS_PetraVector.H"
#include "NLS_PetraGroupInterface.H"
#include "NLS_PetraSharedJacobian.H"

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

  //! Constructor.
  NLS_PetraGroup(Epetra_Vector& x, NLS_PetraSharedJacobian& sj, NLS_PetraGroupInterface& I);

  //! Copy constructor.
  /*! Deep copy. Takes ownership of shared Jacobian if the Jacobian is valid. */
  NLS_PetraGroup(const NLS_PetraGroup& copyFrom);

  //! Deconstructor
  virtual ~NLS_PetraGroup();

  virtual NLS_Group& operator=(const NLS_Group& copyFrom);

  //! See above
  virtual NLS_Group& operator=(const NLS_PetraGroup& copyFrom);

  /** @name "Compute" functions */
  //@{

  virtual const NLS_Vector& computeX(const NLS_Group& grp, const NLS_Vector& d, double step);

  //! See above
  virtual const NLS_Vector& computeX(const NLS_PetraGroup& grp, const NLS_PetraVector& d, double step);

  virtual const NLS_Vector& computeRHS();

  virtual void computeJacobian();

  virtual const NLS_Vector& computeGrad();

  virtual const NLS_Vector& computeNewton(NLS_ParameterList& params);

  virtual double computeNormPredictedRHS(const NLS_Vector& d, double step = 1.0);  

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

  virtual const NLS_Vector& getX() const;

  virtual const NLS_Vector& getRHS() const;
  
  virtual double getNormRHS() const;

  virtual const NLS_Vector& getGrad() const;

  virtual const NLS_Vector& getNewton() const;

  //@}

  virtual NLS_Group* clone() const;

protected:

  //! resets the isValid flags to false
  void resetIsValid();

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

  //! Reference to shared Jacobian matrix 
  NLS_PetraSharedJacobian& sharedJacobian;

  //! Reference to interface for filling the RHS (function) and Jacobian
  NLS_PetraGroupInterface& interface;

};


#endif

// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_GENERICGROUP_H
#define _NLS_GENERICGROUP_H

#include "NLS_Group.H"

//! Semi-Concrete implementation of the pure virtual NLS_Group.

class NLS_GenericGroup : public NLS_Group {

public:

  //! NLS_Group constructor
  NLS_GenericGroup(NLS_Vector& x, NLS_Vector& r);

  //! NLS_Group deconstructor
  virtual ~NLS_GenericGroup();

  //@{ \name "Compute" functions
  virtual NLS_Vector& computeX(const NLS_Group& grp, const NLS_Vector& d, double step);
  virtual NLS_Vector& computeRHS() = 0;
  virtual void computeJacobian() = 0;
  virtual NLS_Vector& computeGrad() = 0;
  virtual NLS_Vector& computeNewton() = 0;
  virtual NLS_Vector& computeNewton(NLS_ParameterList& p) = 0;
  //@}

  //@{ \name "Is" functions
  virtual bool isRHS() const;
  virtual bool isJacobian() const;
  virtual bool isGrad() const;
  virtual bool isNewton() const;
  //@}

  //@{ \name "Get" functions
  virtual NLS_Vector& getX() const;
  virtual NLS_Vector& getRHS() const;
  virtual NLS_Vector& getGrad() const;
  virtual NLS_Vector& getNewton() const;
  //@}

  //@{ \name Creating new NLS_Groups

  virtual NLS_Group* newcopy() = 0;

  //@}
  

};


#endif

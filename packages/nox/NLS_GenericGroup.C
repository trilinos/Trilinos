// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_GenericGroup.H"

//! Constructor
NLS_GenericGroup::NLS_GenericGroup() 
{

}

//! Copy constructor
NLS_GenericGroup::NLS_GenericGroup(NLS_Group& copyFrom) 
{

}

//! Create a new group where the new solution vector is grp.x() + step * d
NLS_GenericGroup::NLS_GenericGroup(NLS_Group& grp, 
				     NLS_Vector& d, double step) 
{

}
 
//! NLS_Group deconstructor
NLS_GenericGroup::~NLS_GenericGroup() {

}

//! Copy Constructor
NLS_GenericGroup& NLS_GenericGroup::operator=(NLS_GenericGroup& 
							copyFrom) {

}

//! Compute and return solution vector
NLS_Vector& NLS_GenericGroup::computeX(NLS_Group& x, NLS_Vector& d, double step) {
  
}

//! Compute and return RHS
NLS_Vector& NLS_GenericGroup::computeRHS() {
  
}

//! Compute RHS
void NLS_GenericGroup::computeJacobian() {
  
}

//! Compute and return gradient 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_GenericGroup::computeGrad() {
  
}

//! Compute and return Newton direction 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_GenericGroup::computeNewton() {
  
}

//! Compute and return Newton direction, using desired accuracy for nonlinear solve
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_GenericGroup::computeNewton(string& name, NLS_Parameter& parameter) {

}

/** @name Checks to see if various objects have been computed. 
 *
 * Returns true if the corresponding "compute" function has been
 * called since the last update to the solution vector (via
 * instantiation or computeX). */

bool NLS_GenericGroup::isRHS() { }
bool NLS_GenericGroup::isJacobian() { }
bool NLS_GenericGroup::isGrad() { }
bool NLS_GenericGroup::isNewton() { }

//! Return solution vector
NLS_Vector& NLS_GenericGroup::getX() { }

//! Return rhs (throws an error if RHS has not been computed)
NLS_Vector& NLS_GenericGroup::getRHS() { }

//! Return gradient (throws an error if gradient has not been computed)
NLS_Vector& NLS_GenericGroup::getGrad() { }

//! Return Newton direction (throws an error if newton direction has not been computed)
NLS_Vector& NLS_GenericGroup::getNewton() { }

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

#ifndef NOX_STATUS_ABSRESID_H
#define NOX_STATUS_ABSRESID_H

#include "NOX_Status_Test.H"	// base class

namespace NOX {

namespace Status {

//! Convergence test based on the norm of the residual.
/*! The norm is computed by a call to the method NOX::Abstract::Group::getNormRHS().  Typically, this computes the 2-norm although the user is free to implement any norm they require in the getNormRHS() method.  For convergence, the tolerance must be below the value of "tolerance" specified in the constructor.   

The second argument to the constructor, "norm" allows the user to scale the norm by \f$ \sqrt{\frac{1}{N}} \f$ where \f$ N \f$ is the number of unknowns in the problem.  For simulations where a mesh is involved (i.e. PDEs) this is a better norm to use than the usual L-2 norm of the residual.  This modification makes the convergence criteria independent of mesh/problem size (a property that is critical for scaling studies).  For example, if the 2-norm is the norm being used, the following equations are solved based on the "normScalingType" flag:

if "norm" is set to "UnscaledNorm":
  \f[
  AbsResid = \sqrt{ \sum_{i=1}^N f_i^2 } 
  \f]

if "norm" is set to "ScaledNorm":
  \f[
  ScaledAbsResid = \sqrt{ \frac{1}{N}\sum_{i=1}^N f_i^2 } 
  \f]
*/
class AbsResid : public Test {

public:

  //! Type that determines whether to scale the norm (computed via a call to NOX::Abstract::Group::getNormRHS()) by the problem size.  
  enum NormType {UnscaledNorm, ScaledNorm};

  //! Constructor.
  /*! <ul>
     <li> tolerance - required residual norm tolerance for convergence.
     <li> norm - If set to "UnscaledNorm" the norm from the NOX::Abstract::Group::getNormRHS() method is used "as is". If set to "ScaledNorm" the norm from NOX::Abstract::Group::getNormRHS() is scaled by \f$ \sqrt{\frac{1}{N}} \f$ where \f$ N \f$ is the number of unknowns in the problem.
     </ul> 
   */
  AbsResid(double tolerance, NormType norm = UnscaledNorm);

  //! Destructor.
  virtual ~AbsResid();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Tolerance required for convergence.
  double tol;

  //! %Status
  StatusType status;

  //! Flag that determines if the norm will be scaled by the problem size.
  NormType normType;

};

} // namespace Status
} // namespace NOX

#endif

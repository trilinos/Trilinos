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

#ifndef NOX_STATUS_SIZEINDABSRESID_H
#define NOX_STATUS_SIZEINDABSRESID_H

#include "NOX_Status_Test.H"	// base class

namespace NOX {

namespace Status {

//! Convergence test based on the absolute norm of the residual scaled by the problem size.  
/*! For simulations where a mesh is involved (i.e. PDEs) this is a better norm to use than a pure L-2 norm of the residual.  This modification makes the abslolute tolerance independent of mesh size. For convergence, the tolerance must be below the value specified in the constructor.  
  
  \f[
  ScaledAbsResid = \sqrt{ \frac{1}{N}\sum_{i=1}^N f_i^2 } 
  \f]

where \f$ f \f$ is the residual vector and \f$ N \f$ is the number of unknowns in the problem.
 */
class SizeIndAbsResid : public Test {

public:

  //! Constructor.
  SizeIndAbsResid(double tolerance);

  //! Destructor.
  virtual ~SizeIndAbsResid();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Tolerance for convergence.
  double tol;

  //! %Status
  StatusType status;

};

} // namespace Status
} // namespace NOX

#endif

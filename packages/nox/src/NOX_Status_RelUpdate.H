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

#ifndef NOX_STATUS_RELUPDATE_H
#define NOX_STATUS_RELUPDATE_H

#include "NOX_Status_Test.H"	// base class

// Forward declarations
namespace NOX {
namespace Abstract {
class Vector;
}
}

namespace NOX {

namespace Status {

//! Convergence test based on the update relative to the solution
/*!
  Convergence test that checks whether or not the following holds:

  \f[
  \frac{1}{N} \sum_{i=1}^N \frac {| x_i - x_{i-1}|}{\epsilon_r |x_i| + \epsilon_a} < 1
  \f]

  Here, \f$\epsilon_r\f$ and \f$\epsilon_a\f$ are specified in the constructor,
  \f$s_i\f$ is the i-th component of the update and \f$x_i\f$ is the
  i-th component of the solution.

*/
class RelUpdate : public Test {

public:

  //! Constructor.
  RelUpdate(double epsilon_r, double epsilon_a);

  //! Destructor.
  virtual ~RelUpdate();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Tolerance for convergence.
  double epsilon_r;

  //! Tolerance for convergence.
  double epsilon_a;

  //! %Status
  StatusType status;

  //! Temporary vector used in computation
  Abstract::Vector* u;

  //! Temporary vector used in computation
  Abstract::Vector* v;

};

} // namespace Status
} // namespace NOX

#endif

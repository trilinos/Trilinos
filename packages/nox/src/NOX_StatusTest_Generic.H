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

#ifndef NOX_STATUSTEST_GENERIC_H
#define NOX_STATUSTEST_GENERIC_H

#include "NOX_Common.H"		// for ostream

namespace NOX {

namespace Solver
{
class Generic;			// forward declaration (included below)
}

//! %Status checkers.
namespace StatusTest {

  //! Status type.
  enum StatusType { Unconverged = 0, Converged = 1, Failed = -1 };

//! Generic status test.
class Generic {

public:

  //! Constructor.
  Generic() {};

  //! Destructor.
  virtual ~Generic() {};

  //! %Test stopping criterion, given a generic problem.
  virtual StatusType checkStatus(const Solver::Generic& problem) = 0;

  //! Return the result of the most recent checkStatus call
  virtual StatusType getStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual ostream& print(ostream& stream, int indent = 0) const = 0;

};

} // namespace StatusTest
} // namespace NOX

ostream& operator<<(ostream& os, NOX::StatusTest::StatusType type);


#endif


// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_STATUS_TEST_H
#define _NOX_STATUS_TEST_H

#include <iostream>

namespace NOX {

namespace Solver
{
class Generic;			// forward declaration
}

//! Status checkers.
namespace Status {

  //! Status type.
  enum StatusType { Unconverged = 0, Converged = 1, Failed = -1};

//! Generic status test.
class Test {

public:

  //! Constructor.
  Test() {};

  //! Destructor.
  virtual ~Test() {};

  //! Test stopping criterion, given a generic problem.
  virtual StatusType operator()(const Solver::Generic& problem) const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual ostream& print(ostream& stream, int indent = 0) const = 0;

};

}
}
#endif

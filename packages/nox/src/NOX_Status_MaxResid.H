// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_STATUS_MAXRESID_H
#define _NOX_STATUS_MAXRESID_H

#include "NOX_Status_Test.H"	// base class

namespace NOX {

namespace Status {

//! Convergence test based on the max-norm of the RHS.
class MaxResid : public Test {

public:

  //! Constructor.
  MaxResid(double tolerance);

  //! Destructor.
  virtual ~MaxResid();

  //! Test stopping criterion, given a generic problem.
  virtual StatusType operator()(const Solver::Generic& problem) const;

  //! Output description of stopping test to output stream
  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  double tol;

};

}
}

#endif

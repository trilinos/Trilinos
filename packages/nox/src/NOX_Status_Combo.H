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

#ifndef NOX_STATUS_COMBO_H
#define NOX_STATUS_COMBO_H

#include "NOX_Status_Test.H"	// base class

namespace NOX {
namespace Status {

  //! Combination of tests.
  /*!
    We always runs through all tests, even if we don't need to. This
    is useful so that the user knows which tests have and have not be
    satisfied.

    In the AND combination, the result is UNCONVERGED if any of the
    tests is UNCONVERGED. Otherwise, the result is equal to the result
    of the first test in the list (either CONVERGED or FAILED). In
    general, one does not imagine mixing CONVERGED and FAILED tests in
    the same AND combination.

    In the OR combination, the result is UNCONVERGED if all of the
    tests are UNCONVERGED. Otherwise, it is the result of the first
    test in the list that is either CONVERGED or FAILED. Therefore, it
    will generally make sense to put the FAILED-type tests at the end
    of the OR list.
  */
class Combo : public Test {

public:

  /*! \brief The test can be either the AND of all the component
    tests, or the OR of all the component tests.
  */
  enum ComboType {AND, OR};

  //! Constructor with a single test.
  Combo(Test& a, ComboType t = OR);

  //! Add another test to this combination.
  virtual Combo& addTest(Test& a);

  //! Destructor
  virtual ~Combo();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;
  
protected:

  //! Use this for operator() when this is an OR type combo
  virtual void orOp(const Solver::Generic& problem);

  //! Use this for operator() when this is an AND type combo
  virtual void andOp(const Solver::Generic& problem);

  //! Check for repeaters.
  bool isSafe(Test& a);

  //! Type of test
  ComboType type;

  //! Vector of tests
  vector<Test*> tests;

  //! %Status
  StatusType status;

};

}
}


#endif

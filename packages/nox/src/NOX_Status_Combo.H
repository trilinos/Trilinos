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
#include "NOX_Common.H"		// class data element (vector)
namespace NOX {
namespace Status {

/*!
  \brief Arbitrary combination of status tests.

  In the \c AND (see NOX::Status::Combo::ComboType) combination, the
  result is \c Unconverged (see NOX::Status::StatusType) if \e any of
  the tests is \c Unconverged. Otherwise, the result is equal to the
  result of the \e first test in the list that is either \c Converged
  or \c Failed. It is not recommended to mix \c Converged and \c
  Failed tests in an \c AND combination.

  In the \c OR combination, the result is \c Unconverged if \e all of
  the tests are \c Unconverged. Otherwise, it is the result of the \e
  first test in the list that is either \c Converged or \c
  Failed. Therefore, it will generally make sense to put the \c Failed
  -type tests at the end of the \c OR list.

  \note We always runs through all tests, even if we don't need
  to. This is useful so that the user knows which tests have and have
  not be satisfied.

  \author Tammy Kolda (SNL 8950)
*/
class Combo : public Test {

public:

  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests.
  */
  enum ComboType {AND, OR};

  //! Constructor
  Combo(ComboType t = OR);

  //! Constructor with a single test.
  Combo(Test& a, ComboType t = OR);

  //! Add another test to this combination.
  virtual Combo& addTest(Test& a);

  //! Destructor
  virtual ~Combo();

  virtual StatusType operator()(const Solver::Generic& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;
  
protected:

  //! Use this for operator() when this is an OR type combo. Updates status.
  virtual void orOp(const Solver::Generic& problem);

  //! Use this for operator() when this is an AND type combo. Updates status.
  virtual void andOp(const Solver::Generic& problem);

  //! Check whether or not it is safe to add a to the list of
  //! tests. This is necessary to avoid any infinite recursions.
  bool isSafe(Test& a);

  //! Type of test
  ComboType type;

  //! Vector of tests
  vector<Test*> tests;

  //! %Status
  StatusType status;

}; // class Combo

} // namespace Status
} // namespace NOX


#endif

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

#ifndef NOX_STATUSTEST_COMBO_H
#define NOX_STATUSTEST_COMBO_H

#include "NOX_StatusTest_Generic.H" // base class
#include "NOX_Common.H"		// class data element (vector)
namespace NOX {
namespace StatusTest {

/*!
  \brief Arbitrary combination of status tests.

  In the \c AND (see NOX::StatusTest::Combo::ComboType) combination, the
  result is \c Unconverged (see NOX::StatusTest::StatusType) if \e any of
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
class Combo : public Generic {

public:

  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests.
  */
  enum ComboType {
    //! Logically "AND" together the results of the tests in this combination
    AND, 
    //! Logically "OR" together the results of the tests in this combination
    OR
  };

  //! Constructor
  Combo(ComboType t);

  //! Constructor with a single test.
  Combo(ComboType t, Generic& a);

  //! Constructor with two tests.
  Combo(ComboType t, Generic& a, Generic& b);

  //! Add another test to this combination. 
  /*!
    Calls isSafe() to determine if it is safe to add \c a to the combination.
  */
  virtual Combo& addStatusTest(Generic& a);

  //! Destructor
  virtual ~Combo();

  /*!
    \brief Tests stopping criterion.

    See addOp() and orOp() for details.
   */
  virtual StatusType checkStatus(const Solver::Generic& problem);

  virtual StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;
  
protected:

  //! Use this for checkStatus when this is an OR type combo. Updates NOX::StatusTest::Combo::status.
  /*!  
    If there is a combination of NOX::StatusTest::Failed and
    NOX::StatusTest::Converged in the tests that are OR'd together,
    the value of status for this test is set to the status of the
    first test it encounters which is not NOX::Status::Unconvered. The
    tests are evaluated in the order that they were added to the
    combination.

    \note We compute the status of all tests in the combination for
    the sake of completeness, even if we could determine the status of
    this combination test without that check.

  */
  virtual void orOp(const Solver::Generic& problem);

  //! Use this for checkStatus when this is an AND type combo. Updates NOX::StatusTest::Combo::status.
  /*!  

    If any tests are NOX::StatusTest::Unconverged, then the status of
    this test is NOX::StatusTest::Unconverged.  If there is a
    combination of NOX::StatusTest::Failed and
    NOX::StatusTest::Converged in the tests that are AND'd together,
    the value of status for this test is set to the status of the
    first test it encounters.  The tests are evaluated in the
    order that they were added to the combination.

    \note We compute the status of all tests in the combination for
    the sake of completeness, even if we could determine the status of
    this combination test without that check.
  */
  virtual void andOp(const Solver::Generic& problem);

  /*! \brief Check whether or not it is safe to add \c a to this list
   of tests.

   This is necessary to avoid any infinite recursions 
   (i.e., a test cannot own a copy of itself).
  */
  bool isSafe(Generic& a);

private:

  //! Type of test
  const ComboType type;

  //! Vector of generic status tests
  vector<Generic*> tests;

  //! %Status
  StatusType status;

}; // class Combo

} // namespace Status
} // namespace NOX


#endif

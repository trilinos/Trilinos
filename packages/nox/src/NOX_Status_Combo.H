// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_STATUS_COMBO_H
#define _NOX_STATUS_COMBO_H

#include "NOX_Status_Test.H"	// base class
#include "NOX_Solver_Generic.H"
#include <vector>

namespace NOX {
namespace Status {

  //! Combination of tests.
class Combo : public Test {

public:

  /*! \brief The test can be either the AND of all the component
   *  tests, or the OR of all the component tests. */
  enum ComboType {AND, OR};

  //! Constructor with a single test.
  Combo(Test& a, ComboType t = OR);

  //! Add another test to this combination.
  virtual Combo& addTest(Test& a);

  //! Destructor
  virtual ~Combo();

  virtual StatusType operator()(const Solver::Generic& problem) const;

  virtual ostream& print(ostream& stream, int indent = 0) const;
  
protected:

  virtual StatusType orOp(const Solver::Generic& problem) const;
  virtual StatusType andOp(const Solver::Generic& problem) const;

  //! Check for repeaters.
  bool isSafe(Test& a);

  //! Type of test
  ComboType type;

  //! Vector of tests
  vector<Test*> tests;

};

}
}


#endif

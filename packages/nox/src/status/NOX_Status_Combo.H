// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_STATUS_BASE_H
#ifndef _NOX_STATUS_BASE_H

#include "NOX_Status_Base.H"	// base class

namespace NOX {
namespace Status {

class Combo : public Test {

public:

  /*! \brief The test can be either the AND of all the component
   *  tests, or the OR of all the component tests. */
  enum ComboType {AND, OR};

  //! Constructor with a single test.
  Combo(Test& a, ComboType t = OR);
  //! Constructor with two tests.
  Combo(Test& a, Test& b, ComboType t = OR);

  //! Add a new test/
  virtual AddTest(Test& a);

  //! Destructor
  virtual ~Combo();


  virtual StopType operator()(const Solver::Base& problem);

  virtual ostream& print(ostream& stream, int indent = 0) const;
  
private:

  //! Type of test
  ComboType type;

  //! Vector of tests
  vector<Test*> tptr;

};

}
}


#endif

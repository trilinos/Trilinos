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

  //! Use this for operator() when this is an OR type combo
  virtual StatusType orOp(const Solver::Generic& problem) const;

  //! Use this for operator() when this is an AND type combo
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

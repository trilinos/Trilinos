// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_StatusTest_Combo.H"
#include "NOX_Utils.H"

NOX::StatusTest::Combo::
Combo(ComboType t, const NOX::Utils* u) :
  type(t)
{
  if (u != NULL)
    utils = *u;

  status = Unevaluated;

}

NOX::StatusTest::Combo::
Combo(ComboType t, 
      const Teuchos::RCP<Generic>& a, 
      const NOX::Utils* u) :
  type(t)
{
  if (u != NULL)
    utils = *u;

  tests.push_back(a);
  status = Unevaluated;

}

NOX::StatusTest::Combo::
Combo(ComboType t, 
      const Teuchos::RCP<Generic>& a, 
      const Teuchos::RCP<Generic>& b, 
      const NOX::Utils* u) :
  type(t)
{
  if (u != NULL)
    utils = *u;

  tests.push_back(a);
  this->addStatusTest(b);
  status = Unevaluated;
}

NOX::StatusTest::Combo& NOX::StatusTest::Combo::
addStatusTest(const Teuchos::RCP<Generic>& a)
{
  if (isSafe(*(a.get())))
    tests.push_back(a);
  else 
  {
    const int indent = 2;
    utils.err() << "\n*** WARNING! ***\n";
    utils.err() << "This combo test currently consists of the following:\n";
    this->print(utils.err(), indent);
    utils.err() << "Unable to add the following test:\n";
    a->print(utils.err(), indent);
    utils.err() << "\n";
  }
  return *this;
}

bool NOX::StatusTest::Combo::isSafe(Generic& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (vector<Teuchos::RCP<Generic> >::iterator i = tests.begin(); i != tests.end(); ++i) 
  {
    
    Combo* ptr = dynamic_cast<Combo*>(i->get());
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

NOX::StatusTest::Combo::~Combo()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::Combo::
checkStatus(const Solver::Generic& problem, 
	    NOX::StatusTest::CheckType checkType)
{
  if (type == OR)
    orOp(problem, checkType);
  else
    andOp(problem, checkType);

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::Combo::getStatus() const
{
  return status;
}

void NOX::StatusTest::Combo::orOp(const Solver::Generic& problem, 
				  NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
  {
    NOX::StatusTest::StatusType s = (*i)->checkStatus(problem, checkType);

    if ((status == Unconverged) && (s != Unconverged)) 
    {
      status = s;

      // Turn off checking for the remaining tests
      if (checkType == NOX::StatusTest::Minimal)
	checkType = NOX::StatusTest::None;
    }

  }

  return;
}

void NOX::StatusTest::Combo::andOp(const Solver::Generic& problem, 
				   NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  bool isUnconverged = false;

  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    NOX::StatusTest::StatusType s = (*i)->checkStatus(problem, checkType);

    // If any of the tests are unconverged, then the AND test is
    // unconverged.
    if (s == Unconverged) 
    {
      isUnconverged = true;
      status = Unconverged;

      // Turn off checking for the remaining tests
      if (checkType == NOX::StatusTest::Minimal)
	checkType = NOX::StatusTest::None;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) && (status == Unconverged)) 
    {
      status = s;
    }

  }

  return;
}


ostream& NOX::StatusTest::Combo::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << endl;

  for (vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
    (*i)->print(stream, indent+2);
    
  return stream;
}

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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_StatusTest_SafeCombo.H"
#include "NOX_Utils.H"

NOX::StatusTest::SafeCombo::SafeCombo(ComboType t) :
  type(t)
{
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo::SafeCombo(ComboType t, 
                                      const Teuchos::RefCountPtr<Generic>& a) :
  type(t)
{
  tests.push_back(a);
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo::SafeCombo(ComboType t, 
                                      const Teuchos::RefCountPtr<Generic>& a, 
                                      const Teuchos::RefCountPtr<Generic>& b) :
  type(t)
{
  tests.push_back(a);
  addStatusTest(b);
  status = Unevaluated;
}

NOX::StatusTest::SafeCombo& NOX::StatusTest::SafeCombo::addStatusTest(const Teuchos::RefCountPtr<Generic>& a)
{
  if (isSafe(a))
    tests.push_back(a);
  else 
  {
    const int indent = 2;
    cout << "\n*** WARNING! ***\n";
    cout << "This combo test currently consists of the following:\n";
    this->print(cout, indent);
    cout << "Unable to add the following test:\n";
    a->print(cout, indent);
    cout << "\n";
  }
  return *this;
}

bool NOX::StatusTest::SafeCombo::isSafe(const Teuchos::RefCountPtr<Generic>& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (a.get() == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (vector<Teuchos::RefCountPtr<Generic> >::iterator i = tests.begin(); 
       i != tests.end(); ++i) 
  {
    
    SafeCombo* ptr = dynamic_cast<SafeCombo*>((*i).get());
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

NOX::StatusTest::SafeCombo::~SafeCombo()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo::checkStatus(const Solver::Generic& problem)
{
  return checkStatusEfficiently(problem, NOX::StatusTest::Minimal);
}

NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo::checkStatusEfficiently(const Solver::Generic& problem, 
					       NOX::StatusTest::CheckType checkType)
{
  if (type == OR)
    orOp(problem, checkType);
  else
    andOp(problem, checkType);

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::SafeCombo::getStatus() const
{
  return status;
}

void NOX::StatusTest::SafeCombo::orOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (vector<Teuchos::RefCountPtr<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
  {
    NOX::StatusTest::StatusType s = (*i)->checkStatusEfficiently(problem, checkType);

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

void NOX::StatusTest::SafeCombo::andOp(const Solver::Generic& problem, NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
    status = Unevaluated;
  else
    status = Unconverged;

  bool isUnconverged = false;

  for (vector<Teuchos::RefCountPtr<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    NOX::StatusTest::StatusType s = (*i)->checkStatusEfficiently(problem, checkType);

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


ostream& NOX::StatusTest::SafeCombo::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
//   stream << setiosflags(ios::left) << setw(13) << setfill('.');
//   if (status == Unconverged) 
//     stream << "**";
//   else if (status == Failed)
//     stream << "Failed";
//   else
//     stream << "Converged";
  stream << status;
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << endl;

  for (vector<Teuchos::RefCountPtr<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) 
    (*i)->print(stream, indent+2);
    
  return stream;
}

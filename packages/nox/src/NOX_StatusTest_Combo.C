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

#include "NOX_StatusTest_Combo.H"
#include "NOX_Utils.H"

using namespace NOX::StatusTest;

Combo::Combo(ComboType t) :
  type(t)
{
  status = Unconverged;
}

Combo::Combo(ComboType t, Generic& a) :
  type(t)
{
  tests.push_back(&a);
  status = Unconverged;
}

Combo::Combo(ComboType t, Generic& a, Generic& b) :
  type(t)
{
  tests.push_back(&a);
  addStatusTest(b);
  status = Unconverged;
}

Combo& Combo::addStatusTest(Generic& a)
{
  if (isSafe(a))
    tests.push_back(&a);
  else 
    if (Utils::doPrint(Utils::Warning)) {
      const int indent = 2;
      cout << "\n*** WARNING! ***\n";
      cout << "This combo test currently consists of the following:\n";
      this->print(cout, indent);
      cout << "Unable to add the following test:\n";
      a.print(cout, indent);
      cout << "\n";
    }
  return *this;
}

bool Combo::isSafe(Generic& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (vector<Generic*>::iterator i = tests.begin(); i != tests.end(); ++i) {
    
    Combo* ptr = dynamic_cast<Combo*>(*i);
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

Combo::~Combo()
{
}

StatusType Combo::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;

  if (type == OR)
    orOp(problem);
  else
    andOp(problem);

  return status;
}

StatusType Combo::getStatus() const
{
  return status;
}

void Combo::orOp(const Solver::Generic& problem)
{
  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (vector<Generic*>::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    StatusType s = (*i)->checkStatus(problem);

    if ((status == Unconverged) && (s != Unconverged)) {
      status = s;
    }

  }

  return;
}

void Combo::andOp(const Solver::Generic& problem)
{
  bool isUnconverged = false;

  for (vector<Generic*>::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    StatusType s = (*i)->checkStatus(problem);

    // If any of the tests are unconverged, then the AND test is
    // unconverged.
    if (s == Unconverged) {
      isUnconverged = true;
      status = Unconverged;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) && (status == Unconverged)) {
      status = s;
    }

  }

  return;
}


ostream& Combo::print(ostream& stream, int indent) const
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

  for (vector<Generic*>::const_iterator i = tests.begin(); i != tests.end(); ++i) 
    (*i)->print(stream, indent+2);
    
  return stream;
}

// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Status_Combo.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX::Status;

Combo::Combo(Test& a, ComboType t)
{
  type = t;
  tests.push_back(&a);
}

Combo& Combo::addTest(Test& a)
{
  if (isSafe(a))
    tests.push_back(&a);
  else 
    if (Utils::doPrint(1)) {
      cout << "\n*** WARNING! ***\n";
      cout << "This combo tests currently consists of the following:\n";
      this->print(cout, 2);
      cout << "Unable to add the following test:\n";
      a.print(cout, 2);
      cout << "\n";
    }
  return *this;
}

bool Combo::isSafe(Test& a)
{
  if (&a == this)
    return false;
  
  Combo* ptr;
  for (vector<Test*>::iterator i = tests.begin(); i != tests.end(); ++i) {
    
    ptr = dynamic_cast<Combo*>(*i);
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  return true;
}

Combo::~Combo()
{
}

StatusType Combo::operator()(const Solver::Generic& problem) const
{
  if (type == OR)
    return orOp(problem);
  else
    return andOp(problem);
}

StatusType Combo::orOp(const Solver::Generic& problem) const
{

  for (vector<Test*>::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    StatusType status = (*i)->operator()(problem);
    if (status != Unconverged) 
      return status;

  }

  return Unconverged;
}

StatusType Combo::andOp(const Solver::Generic& problem) const
{
  vector<Test*>::const_iterator i = tests.begin();
  StatusType status = (*i)->operator()(problem);

  if (status == Unconverged)
    return Unconverged;

  StatusType newstatus;

  for (; i != tests.end(); ++i) 
    if (status != (*i)->operator()(problem))
      return Unconverged;

  return status;
}


ostream& Combo::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "Combo ";
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Test:" << endl;

  for (vector<Test*>::const_iterator i = tests.begin(); i != tests.end(); ++i) 
    (*i)->print(stream, indent+2);
    
  return stream;
}

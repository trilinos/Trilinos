// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  for (std::vector<Teuchos::RCP<Generic> >::iterator i = tests.begin(); i != tests.end(); ++i)
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
  for (std::vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i)
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

  for (std::vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i) {

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


std::ostream& NOX::StatusTest::Combo::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << std::endl;

  for (std::vector<Teuchos::RCP<Generic> >::const_iterator i = tests.begin(); i != tests.end(); ++i)
    (*i)->print(stream, indent+2);

  return stream;
}

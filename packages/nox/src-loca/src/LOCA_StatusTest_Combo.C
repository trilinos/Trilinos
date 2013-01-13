// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "LOCA_StatusTest_Combo.H"

// FIXME Really necessary?
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

LOCA::StatusTest::Combo::
Combo(ComboType t, const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr ) :
  type(t)
{
  if ( globalDataPtr.is_valid_ptr() && !globalDataPtr.is_null() )
    globalDataPtr_ = globalDataPtr;

  status = LOCA::StatusTest::Unevaluated;
}

LOCA::StatusTest::Combo::
Combo(ComboType t, 
      const Teuchos::RCP<Abstract>& a,
      const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr ) :
  type(t)
{
  if ( globalDataPtr.is_valid_ptr() && !globalDataPtr.is_null() )
    globalDataPtr_ = globalDataPtr;

  tests.push_back(a);

  status = LOCA::StatusTest::Unevaluated;
}

LOCA::StatusTest::Combo::
Combo(ComboType t, 
      const Teuchos::RCP<Abstract>& a,
      const Teuchos::RCP<Abstract>& b,
      const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr ) :
  type(t)
{
  if ( globalDataPtr.is_valid_ptr() && !globalDataPtr.is_null() )
    globalDataPtr_ = globalDataPtr;

  tests.push_back(a);
  this->addStatusTest(b);

  status = LOCA::StatusTest::Unevaluated;
}

LOCA::StatusTest::Combo& LOCA::StatusTest::Combo::
addStatusTest(const Teuchos::RCP<Abstract>& a)
{
  if (isSafe(*(a.get())))
    tests.push_back(a);
  else 
  {
    const int indent = 2;
    globalDataPtr_->locaUtils->err() << "\n*** WARNING! ***\n";
    globalDataPtr_->locaUtils->err() << "This combo test currently consists of the following:\n";
    this->print(globalDataPtr_->locaUtils->err(), indent);
    globalDataPtr_->locaUtils->err() << "Unable to add the following test:\n";
    a->print(globalDataPtr_->locaUtils->err(), indent);
    globalDataPtr_->locaUtils->err() << "\n";
  }
  return *this;
}

bool LOCA::StatusTest::Combo::isSafe(Abstract& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (std::vector<Teuchos::RCP<Abstract> >::iterator i = tests.begin(); i != tests.end(); ++i)
  {
    Combo* ptr = dynamic_cast<Combo*>(i->get());
    if (ptr != NULL)
      if (!ptr->isSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

LOCA::StatusTest::Combo::~Combo()
{
}

LOCA::StatusTest::StatusType LOCA::StatusTest::Combo::
//checkStatus(const LOCA::Stepper& stepper,
checkStatus(const LOCA::Abstract::Iterator& stepper,
	    LOCA::StatusTest::CheckType checkType)
{
  if (type == OR)
    orOp(stepper, checkType);
  else
    andOp(stepper, checkType);

  return status;
}

LOCA::StatusTest::StatusType LOCA::StatusTest::Combo::
getStatus() const
{
  return status;
}

//void LOCA::StatusTest::Combo::orOp(const LOCA::Stepper& stepper,
void LOCA::StatusTest::Combo::orOp(const LOCA::Abstract::Iterator& stepper,
				   LOCA::StatusTest::CheckType checkType)
{
  if (checkType == LOCA::StatusTest::None)
    status = LOCA::StatusTest::Unevaluated;
  else
    status = LOCA::StatusTest::NotFinished;

  // Checks the status of each test. The first test it encounters, if
  // any, that is NotFinished is the status that it sets itself too.
  for (std::vector<Teuchos::RCP<Abstract> >::const_iterator i = tests.begin(); i != tests.end(); ++i)
  {
    LOCA::StatusTest::StatusType s = (*i)->checkStatus(stepper, checkType);

    if ((status == LOCA::StatusTest::NotFinished) && (s != LOCA::StatusTest::NotFinished))
    {
      status = s;

      // Turn off checking for the remaining tests
      if (checkType == LOCA::StatusTest::Minimal)
	checkType = LOCA::StatusTest::None;
    }

  }

  return;
}

//void LOCA::StatusTest::Combo::andOp(const LOCA::Stepper& stepper,
void LOCA::StatusTest::Combo::andOp(const LOCA::Abstract::Iterator& stepper,
				   LOCA::StatusTest::CheckType checkType)
{
  if (checkType == LOCA::StatusTest::None)
    status = LOCA::StatusTest::Unevaluated;
  else
    status = LOCA::StatusTest::NotFinished;

  bool isUnconverged = false;

  for (std::vector<Teuchos::RCP<Abstract> >::const_iterator i = tests.begin(); i != tests.end(); ++i) {

    LOCA::StatusTest::StatusType s = (*i)->checkStatus(stepper, checkType);

    // If any of the tests are NotFinished, then the AND test is
    // NotFinished.
    if (s == LOCA::StatusTest::NotFinished)
    {
      isUnconverged = true;
      status = LOCA::StatusTest::NotFinished;

      // Turn off checking for the remaining tests
      if (checkType == LOCA::StatusTest::Minimal)
	checkType = LOCA::StatusTest::None;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) && (status == LOCA::StatusTest::NotFinished))
    {
      status = s;
    }

  }

  return;
}


std::ostream& LOCA::StatusTest::Combo::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << ((type == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << std::endl;

  for (std::vector<Teuchos::RCP<Abstract> >::const_iterator i = tests.begin(); i != tests.end(); ++i)
    (*i)->print(stream, indent+2);
    
  return stream;
}


//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "AztecOO_StatusTestCombo.h"

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t)
  : AztecOO_StatusTest() {
  type_ = t;
  status_ = Unchecked;
}

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a)
  : AztecOO_StatusTest() {
  type_ = t;
  tests_.push_back(&a);
  status_ = Unchecked;
}

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b)
  : AztecOO_StatusTest() {
  type_ = t;
  tests_.push_back(&a);
  AddStatusTest(b);
  status_ = Unchecked;
}

AztecOO_StatusTestCombo& AztecOO_StatusTestCombo::AddStatusTest(AztecOO_StatusTest& a)
{
  if (IsSafe(a))
    tests_.push_back(&a);
  else 
    // if (Utils::doPrint(Utils::Warning)) 
    {
      const int indent = 2;
      cout << "\n*** WARNING! ***\n";
      cout << "This combo test currently consists of the following:\n";
      this->Print(cout, indent);
      cout << "Unable to add the following test:\n";
      a.Print(cout, indent);
      cout << "\n";
    }
  return *this;
}

bool AztecOO_StatusTestCombo::IsSafe(AztecOO_StatusTest& a)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&a == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (std::vector<AztecOO_StatusTest*>::iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    AztecOO_StatusTestCombo* ptr = dynamic_cast<AztecOO_StatusTestCombo*>(*i);
    if (ptr != NULL)
      if (!ptr->IsSafe(a))
	return false;
  }

  // Otherwise, it's safe to add a to the list.
  return true;
}

bool AztecOO_StatusTestCombo::ResidualVectorRequired() const
{
  // If any of the StatusTest object require the residual vector, then return true.

  // Recursively test this property.
  for (std::vector<AztecOO_StatusTest * >::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    AztecOO_StatusTest* ptr = dynamic_cast<AztecOO_StatusTest*>(*i);
    if (ptr != NULL)
      if (ptr->ResidualVectorRequired())
	return true;
  }

  // Otherwise we don't need residual vector.
  return false;
}

AztecOO_StatusType AztecOO_StatusTestCombo::CheckStatus(int CurrentIter, 
							      Epetra_MultiVector * CurrentResVector, 
							      double CurrentResNormEst,  
							      bool SolutionUpdated) {
  status_ = Unconverged;

  if (type_ == OR)
    OrOp(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);
  else if (type_ == AND)
    AndOp(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);
  else
    SeqOp(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

  return status_;
}


void AztecOO_StatusTestCombo::OrOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, 
				   double CurrentResNormEst, bool SolutionUpdated) {

  bool isFailed = false;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (std::vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    AztecOO_StatusType s = (*i)->CheckStatus(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);
    
    // Check for failure and NaN.  Combo treats NaNs as Fails
    if (s==Failed || s==NaN) isFailed = true;

    if ((status_ == Unconverged) && (s != Unconverged)) {
      status_ = s;
    }

  }

    // Any failure is a complete failure
    if (isFailed) status_ = Failed;

  return;
}

void AztecOO_StatusTestCombo::AndOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
				    bool SolutionUpdated) {
  bool isUnconverged = false;
  bool isFailed = false;

  for (std::vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    AztecOO_StatusType s = (*i)->CheckStatus(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

    // Check for failure and NaN.  Combo treats NaNs as Fails
    if (s==Failed || s==NaN) isFailed = true;

    // If any of the tests are unconverged, then the AND test is
    // unconverged.
    if (s == Unconverged) {
      isUnconverged = true;
      status_ = Unconverged;
    }

    // If this is the first test and it's converged/failed, copy its
    // status to the combo status.
    if ((!isUnconverged) && (status_ == Unconverged)) {
      status_ = s;
    }

  }

    // Any failure is a complete failure
    if (isFailed) status_ = Failed;

  return;
}

void AztecOO_StatusTestCombo::SeqOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
				    bool SolutionUpdated) {

  for (std::vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    AztecOO_StatusType s = (*i)->CheckStatus(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

    // Check for failure and NaN.  Combo treats NaNs as Fails
    if (s==Failed || s==NaN) {
      status_ = Failed;
      return;
    }
    else if (s==Unconverged) {
      status_ = s;
      return;
    }
  }

  // If we make it here, we have converged
  status_ = Converged;

  return;
}


ostream& AztecOO_StatusTestCombo::Print(ostream& stream, int indent) const {
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  PrintStatus(stream, status_);
  stream << ((type_ == OR) ? "OR" : (type_ == AND) ? "AND" :"SEQ");
  stream << " Combination";
  stream << " -> " << endl;

  for (std::vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    (*i)->Print(stream, indent+2);
    
  return stream;
}


/* Copyright (2001) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 *
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "AztecOO_StatusTestCombo.h"

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t)
  : AztecOO_StatusTest() {
  type_ = t;
  status_ = Unconverged;
}

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a)
  : AztecOO_StatusTest() {
  type_ = t;
  tests_.push_back(&a);
  status_ = Unconverged;
}

AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t, AztecOO_StatusTest& a, AztecOO_StatusTest& b)
  : AztecOO_StatusTest() {
  type_ = t;
  tests_.push_back(&a);
  AddStatusTest(b);
  status_ = Unconverged;
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
  for (vector<AztecOO_StatusTest*>::iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
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
  for (vector<AztecOO_StatusTest * const>::iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
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
							      bool SolutionUpdated)
{
  status_ = Unconverged;

  if (type_ == OR)
    OrOp(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);
  else
    AndOp(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

  return status_;
}

AztecOO_StatusType AztecOO_StatusTestCombo::GetStatus() const
{
  return status_;
}

void AztecOO_StatusTestCombo::OrOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
			   bool SolutionUpdated)
{
  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    AztecOO_StatusType s = (*i)->CheckStatus(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

    if ((status_ == Unconverged) && (s != Unconverged)) {
      status_ = s;
    }

  }

  return;
}

void AztecOO_StatusTestCombo::AndOp(int CurrentIter, Epetra_MultiVector * CurrentResVector, double CurrentResNormEst,
			   bool SolutionUpdated)
{
  bool isUnconverged = false;

  for (vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    AztecOO_StatusType s = (*i)->CheckStatus(CurrentIter, CurrentResVector, CurrentResNormEst, SolutionUpdated);

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

  return;
}


ostream& AztecOO_StatusTestCombo::Print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  PrintStatus(stream, status_);
  stream << ((type_ == OR) ? "OR" : "AND");
  stream << " Combination";
  stream << " -> " << endl;

  for (vector<AztecOO_StatusTest*>::const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    (*i)->Print(stream, indent+2);
    
  return stream;
}

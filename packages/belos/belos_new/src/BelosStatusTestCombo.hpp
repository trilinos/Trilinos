
// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER
//

#ifndef BELOS_STATUS_TEST_COMBO_H
#define BELOS_STATUS_TEST_COMBO_H

/*!
  \file BelosStatusTestCombo.hpp
  \brief Belos::StatusTest for logically combining several status tests.
*/
  
#include "BelosStatusTest.hpp"
#include <vector>

/*! 
  \class Belos::StatusTestCombo
  \brief A class for extending the status testing capabilities of Belos via logical combinations.

  StatusTestCombo is an interface that can be implemented to extend the convergence testing
  capabilities of Belos.  This class supports composite tests.  In this situation,
  two or more existing StatusTestCombo objects test1 and test2 can be used to create a new test.
  For all combinations, if any tests returns Failed or returns not-a-number (NaN) status, then the combination test 
  returns Failed.
  There are three possible combinations:
  <ol>
  <li> OR combination:
  If an OR combination is selected, the status returns Converged if any one of the subtest returns
  as Converged.  
  <li> AND combination:
  If an AND combination is selected, the status returns Converged only when all subtests return as Converged.
  <li> SEQ combination:
  SEQ is a form of AND that will perform subtests in sequence.  If the first test returns Unconverged, Failed or NaN,
  no other subtests are done, and the status is returned as Unconverged if the first test was Unconverged, or as
  Failed if the first test was Failed or NaN.  If the first test returns Converged, the second test is checked in 
  the same fashion as the first.  If the second test is Converged, the third one is tested, and so on.
  
  The purpose of the SEQ combination is to allow the addition of expensive but more rigorous convergence tests.  For
  example, we could define a test that used the implicit residual vector (the one produced by the iterative method)
  as the first subtest and define a second test using the explicitly computed residual vector.  Explicitly computing
  the residual requires a matrix multiplication with the original matrix operator, an expensive operation.  By using
  the SEQ combination, we can avoid the matrix multiplication associated with the explicit residual calculation
  until the implicit residual is small.
  </ol>
*/

namespace Belos {

template <class TYPE, class OP, class MV>
class StatusTestCombo: public StatusTest<TYPE,OP,MV> {
	
 public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef std::vector< StatusTest<TYPE,OP,MV>* > st_vector;
  typedef typename st_vector::iterator iterator;
  typedef typename st_vector::const_iterator const_iterator;

#endif // DOXYGEN_SHOULD_SKIP_THIS
  
  //@{ \name Enums.
  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests, or a sequential AND (SEQ).
  */
  enum ComboType {AND,  /*!< Require all subtests to be satisfied. */
		  OR,   /*!< Require one or the other subtests to be satisfied. */
		  SEQ   /*!< Requires all subtests to be satisfied, but stops check after the first failed 
			  or unconverged status. */
  };
  //@}

  //@{ \name Constructors / Destructor.

  //! Constructor
  StatusTestCombo(ComboType t);

  //! Single test constructor.
  StatusTestCombo(ComboType t, StatusTest<TYPE,OP,MV>& test1);

  //! Dual test constructor.
  StatusTestCombo(ComboType t, StatusTest<TYPE,OP,MV>& test1, StatusTest<TYPE,OP,MV>& test2);

  //! Add another test to this combination.
  StatusTestCombo<TYPE,OP,MV>& AddStatusTest(StatusTest<TYPE,OP,MV>& add_test);

  //! Destructor
  virtual ~StatusTestCombo() {};
  //@}

  //@{ \name Status methods
  
  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met using the current information from the 
    iterative solver.
  */
  StatusType CheckStatus( IterativeSolver<TYPE,OP,MV>* iSolver );

  //! Return the result of the most recent CheckStatus call.
  StatusType GetStatus() const { return(status_); };

  //@}

  //@{ \name Reset methods

  //! Resets all the status tests in this combination to their initial internal state.
  /*! This should be done when the status test is being reused with another solver or linear problem.
  */
  void Reset(); 

  //@}

  //@{ \name Accessor methods

  //! Returns the maximum number of iterations set in the constructor.
  ComboType GetComboType() const {return(type_);};

  //@}

  //@{ \name Attribute methods

  //! Indicates if residual vector is required by this convergence test.
  /*! If this method returns true, then one or more of the StatusTest objects that make up this combined
    test requires the norm of the true residual vector to perform its test.
  */
  bool ResidualVectorRequired() const;

  //@}
  //@{ \name Print methods
  
  //! Output formatted description of stopping test to output stream
  ostream& Print(ostream& os, int indent = 0) const;
  
  //@}

protected:

  //@{ \name Internal methods.
  //! Use this for checkStatus when this is an OR type combo. Updates status.
  void OrOp( IterativeSolver<TYPE,OP,MV>* iSolver );

  //! Use this for checkStatus when this is an AND type combo. Updates status.
  void AndOp( IterativeSolver<TYPE,OP,MV>* iSolver );

  //! Use this for checkStatus when this is a sequential AND type combo. Updates status.
  void SeqOp( IterativeSolver<TYPE,OP,MV>* iSolver );

  //! Check whether or not it is safe to add a to the list of
  //! tests. This is necessary to avoid any infinite recursions.
  bool IsSafe(StatusTest<TYPE,OP,MV>& test1);
  //@}

 private:

  //@{ \name Private data members.
  //! Type of test
  ComboType type_;

  //! Vector of generic status tests
  st_vector tests_;

  //! Status
   StatusType status_;
  //@}

};

template <class TYPE, class OP, class MV>
StatusTestCombo<TYPE,OP,MV>::StatusTestCombo(ComboType t)
{
  type_ = t;
  status_ = Unchecked;
}

template <class TYPE, class OP, class MV>
StatusTestCombo<TYPE,OP,MV>::StatusTestCombo(ComboType t, StatusTest<TYPE,OP,MV>& test1)
{
  type_ = t;
  tests_.push_back(&test1);
  status_ = Unchecked;
}

template <class TYPE, class OP, class MV>
StatusTestCombo<TYPE,OP,MV>::StatusTestCombo(ComboType t, StatusTest<TYPE,OP,MV>& test1, StatusTest<TYPE,OP,MV>& test2)
{
  type_ = t;
  tests_.push_back(&test1);
  AddStatusTest(test2);
  status_ = Unchecked;
}

template <class TYPE, class OP, class MV>
StatusTestCombo<TYPE,OP,MV>& StatusTestCombo<TYPE,OP,MV>::AddStatusTest(StatusTest<TYPE,OP,MV>& add_test)
{
  if (IsSafe(add_test))
    tests_.push_back(&add_test);
  else
    {
      const int indent = 2;
      cout << "\n*** WARNING! ***\n";
      cout << "This combo test currently consists of the following:\n";
      this->Print(cout, indent);
      cout << "Unable to add the following test:\n";
      add_test.Print(cout, indent);
      cout << "\n";
    }
  return *this;
}

template <class TYPE, class OP, class MV>
bool StatusTestCombo<TYPE,OP,MV>::IsSafe(StatusTest<TYPE,OP,MV>& test1)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (&test1 == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    StatusTestCombo<TYPE,OP,MV>* ptr = dynamic_cast< StatusTestCombo<TYPE,OP,MV> *>(*i);
    if (ptr != NULL)
      if (!ptr->IsSafe(test1))
        return false;
  }
  return true;
}

template <class TYPE, class OP, class MV>
bool StatusTestCombo<TYPE,OP,MV>::ResidualVectorRequired() const
{
  // If any of the StatusTest object require the residual vector, then return true.
  
  // Recursively test this property.
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    StatusTest<TYPE,OP,MV>* ptr = dynamic_cast< StatusTest<TYPE,OP,MV> *>(*i);
    if (ptr != NULL)
      if (ptr->ResidualVectorRequired())
        return true;
  }
  
  // Otherwise we don't need residual vector.
  return false;
}

template <class TYPE, class OP, class MV>
StatusType StatusTestCombo<TYPE,OP,MV>::CheckStatus( IterativeSolver<TYPE,OP,MV>* iSolver )
{
  status_ = Unconverged;

  if (type_ == OR)
    OrOp( iSolver );
  else if (type_ == AND)
    AndOp( iSolver );
  else
    SeqOp( iSolver );

  return status_;
}

template <class TYPE, class OP, class MV>
void StatusTestCombo<TYPE,OP,MV>::Reset( )
{
  // Resets all status tests.
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    {
      (*i)->Reset();
    }
  return;
}

template <class TYPE, class OP, class MV>
void StatusTestCombo<TYPE,OP,MV>::OrOp( IterativeSolver<TYPE,OP,MV>* iSolver )
{
  bool isFailed = false;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    {
      StatusType s = (*i)->CheckStatus( iSolver );
      
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

template <class TYPE, class OP, class MV>
void StatusTestCombo<TYPE,OP,MV>::AndOp( IterativeSolver<TYPE,OP,MV>* iSolver )
{
  bool isUnconverged = false;
  bool isFailed = false;
  
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    StatusType s = (*i)->CheckStatus( iSolver );

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

template <class TYPE, class OP, class MV>
void StatusTestCombo<TYPE,OP,MV>::SeqOp( IterativeSolver<TYPE,OP,MV>* iSolver ) 
{
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    StatusType s = (*i)->CheckStatus( iSolver );

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

template <class TYPE, class OP, class MV>
ostream& StatusTestCombo<TYPE,OP,MV>::Print(ostream& os, int indent) const {
  for (int j = 0; j < indent; j ++)
    os << ' ';
  PrintStatus(os, status_);
  os << ((type_ == OR) ? "OR" : (type_ == AND) ? "AND" :"SEQ");
  os << " Combination";
  os << " -> " << endl;

  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i)
    (*i)->Print(os, indent+2);

  return os;
}

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_COMBO_H */

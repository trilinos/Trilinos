
// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  If an OR combination is selected, the status returns Converged if any one of the subtest returns.
  All of the tests are executed, whether or not any of them return Converged.
  as Converged.  
  <li> AND combination:
  If an AND combination is selected, the status returns Converged only when all subtests return as Converged.
  All of the tests are executed, whether or not any of them return Failed.
  <li> SEQ combination:
  SEQ is a form of AND that will perform subtests in sequence.  If the first test returns Passed, Failed or Undefined,
  no other subtests are done, and the status is returned as Failed if the first test was Failed, or as
  Failed if the first test was Failed or NaN.  If the first test returns Converged, the second test is checked in 
  the same fashion as the first.  If the second test is Converged, the third one is tested, and so on.
  
  The purpose of the SEQ combination is to allow the addition of expensive but more rigorous convergence tests.  For
  example, we could define a test that used the implicit residual std::vector (the one produced by the iterative method)
  as the first subtest and define a second test using the explicitly computed residual std::vector.  Explicitly computing
  the residual requires a matrix multiplication with the original matrix operator, an expensive operation.  By using
  the SEQ combination, we can avoid the matrix multiplication associated with the explicit residual calculation
  until the implicit residual is small.
  </ol>
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestCombo: public StatusTest<ScalarType,MV,OP> {
	
 public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef std::vector< Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > st_vector;
  typedef typename st_vector::iterator iterator;
  typedef typename st_vector::const_iterator const_iterator;

#endif // DOXYGEN_SHOULD_SKIP_THIS
  
  //! @name Enums
  //@{ 
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

  //! @name Constructors / Destructor
  //@{ 

  //! Constructor
  StatusTestCombo(ComboType t);

  //! Single test constructor.
  StatusTestCombo(ComboType t, 
		  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test1);

  //! Dual test constructor.
  StatusTestCombo(ComboType t, 
		  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test1, 
		  const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test2);

  /// \brief Add another test to this combination.
  ///
  /// Only add the test if doing so would not in infinite recursion.
  StatusTestCombo<ScalarType,MV,OP>& addStatusTest(const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& add_test);

  //! Destructor
  virtual ~StatusTestCombo() {};
  //@}

  //! @name Status methods
  //@{ 
  
  /// \brief Check convergence status of the iterative solver.
  ///
  /// Return one of the following values: Passed (the convergence
  /// criteria are met), Failed (they are not met) or Undefined (we
  /// can't tell).
  StatusType checkStatus( Iteration<ScalarType,MV,OP>* iSolver );

  /// \brief Return the result of the most recent checkStatus call.
  ///
  /// If checkStatus has not yet been called, return the default status.
  StatusType getStatus() const { return(status_); };

  //@}

  //! @name Reset methods
  //@{ 

  //! Resets all the status tests in this combination to their initial internal state.
  /*! This should be done when the status test is being reused with another solver or linear problem.
  */
  void reset(); 

  //@}

  //! @name Accessor methods
  //@{ 

  //! Return the type of combination (OR, AND, or SEQ).
  ComboType getComboType() const { return type_; }

  //! Return the vector of status tests
  st_vector getStatusTests() { return tests_; }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream
  void print(std::ostream& os, int indent = 0) const;
  
  //@}

protected:

  //! @name Internal methods.
  //@{ 
  //! Use this for checkStatus when this is an OR type combo. Updates status.
  void orOp( Iteration<ScalarType,MV,OP>* iSolver );

  //! Use this for checkStatus when this is an AND type combo. Updates status.
  void andOp( Iteration<ScalarType,MV,OP>* iSolver );

  //! Use this for checkStatus when this is a sequential AND type combo. Updates status.
  void seqOp( Iteration<ScalarType,MV,OP>* iSolver );

  //! Check whether or not it is safe to add a to the list of
  //! tests. This is necessary to avoid any infinite recursions.
  bool isSafe( const Teuchos:: RCP<StatusTest<ScalarType,MV,OP> >& test1);
  //@}

 private:

  //! @name Private data members.
  //@{ 
  //! The type of combination (OR, AND, or SEQ)
  ComboType type_;

  //! Vector of generic status tests
  st_vector tests_;

  /// \brief The current status
  StatusType status_;
  //@}

};

template <class ScalarType, class MV, class OP>
StatusTestCombo<ScalarType,MV,OP>::StatusTestCombo(ComboType t)
{
  type_ = t;
  status_ = Undefined;
}

template <class ScalarType, class MV, class OP>
StatusTestCombo<ScalarType,MV,OP>::StatusTestCombo(ComboType t, 
						   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test1)
{
  type_ = t;
  tests_.push_back(test1);
  status_ = Undefined;
}

template <class ScalarType, class MV, class OP>
StatusTestCombo<ScalarType,MV,OP>::StatusTestCombo(ComboType t, 
						   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test1, 
						   const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test2)
{
  type_ = t;
  tests_.push_back(test1);
  addStatusTest(test2);
  status_ = Undefined;
}

template <class ScalarType, class MV, class OP>
StatusTestCombo<ScalarType,MV,OP>& StatusTestCombo<ScalarType,MV,OP>::addStatusTest(const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& add_test)
{
  if (isSafe(add_test))
    tests_.push_back(add_test);
  else
    {
      const int indent = 2;
      std::cout << "\n*** WARNING! ***\n";
      std::cout << "This combo test currently consists of the following:\n";
      this->print(std::cout, indent);
      std::cout << "Unable to add the following test:\n";
      add_test->print(std::cout, indent);
      std::cout << "\n";
    }
  return *this;
}

template <class ScalarType, class MV, class OP>
bool StatusTestCombo<ScalarType,MV,OP>::isSafe( const Teuchos::RCP<StatusTest<ScalarType,MV,OP> >& test1)
{
  // Are we trying to add "this" to "this"? This would result in an infinite recursion.
  if (test1.get() == this)
    return false;
  
  // Recursively test that we're not adding something that's already
  // in the list because that can also lead to infinite recursions.
  for (iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    StatusTestCombo<ScalarType,MV,OP>* ptr = dynamic_cast<StatusTestCombo<ScalarType,MV,OP> *>(i->get());
    if (ptr != NULL)
      if (!ptr->isSafe(test1))
        return false;
  }
  return true;
}

template <class ScalarType, class MV, class OP>
StatusType StatusTestCombo<ScalarType,MV,OP>::checkStatus( Iteration<ScalarType,MV,OP>* iSolver )
{
  status_ = Failed;

  if (type_ == OR)
    orOp( iSolver );
  else if (type_ == AND)
    andOp( iSolver );
  else
    seqOp( iSolver );

  return status_;
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::reset( )
{
  // Resets all status tests in my list.
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    {
      (*i)->reset();
    }
  // Reset my status.
  status_ = Undefined;
  //
  return;
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::orOp( Iteration<ScalarType,MV,OP>* iSolver )
{
  status_ = Failed;

  // Checks the status of each test. The first test it encounters, if
  // any, that is unconverged is the status that it sets itself too.
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) 
    {
      StatusType s = (*i)->checkStatus( iSolver );

      // Check for failure.
      if (s==Passed) status_ = Passed;
    }
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::andOp( Iteration<ScalarType,MV,OP>* iSolver )
{
  bool isFailed = false;
  
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) {
    
    StatusType s = (*i)->checkStatus( iSolver );

    // Check for failure.
    if (s==Failed) isFailed = true;

    // If any of the tests are failed, then the AND test is failed.
    if (s == Failed) {
      status_ = Failed;
    }

    // If this is the first test and it's failed, copy its
    // status to the combo status.
    if ((!isFailed) && (status_ == Failed)) {
      status_ = s;
    }
  }
  
  // Any failure is a complete failure
  if (isFailed) status_ = Failed;
  
  return;
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::seqOp( Iteration<ScalarType,MV,OP>* iSolver ) 
{
  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i) {

    StatusType s = (*i)->checkStatus( iSolver );

    // Check for failure.
    if (s==Failed) {
      status_ = Failed;
      return;
    }
    else if (s==Undefined) {
      status_ = s;
      return;
    }
  }
  // If we make it here, we have converged
  status_ = Passed;

  return;
}

template <class ScalarType, class MV, class OP>
void StatusTestCombo<ScalarType,MV,OP>::print(std::ostream& os, int indent) const {
  for (int j = 0; j < indent; j ++)
    os << ' ';
  this->printStatus(os, status_);
  os << ((type_ == OR) ? "OR" : (type_ == AND) ? "AND" :"SEQ");
  os << " Combination";
  os << " -> " << std::endl;

  for (const_iterator i = tests_.begin(); i != tests_.end(); ++i)
    (*i)->print(os, indent+2);
}

} // end namespace Belos

#endif /* BELOS_STATUS_TEST_COMBO_H */

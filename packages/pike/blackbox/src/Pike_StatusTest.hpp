#ifndef PIKE_STATUS_TEST_HPP
#define PIKE_STATUS_TEST_HPP

#include "Pike_BlackBox_config.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <iostream>

namespace pike {

  static const int defaultIndentation = 3;
  static const int statusIndentation = 13;

  class Solver;

  //! The status of a solver
  enum SolveStatus {
    //! The solver is converged.  This triggers termination of the solver.
    CONVERGED,
    //! The solver has failed.  This triggers termination of the solver.
    FAILED,
    //! The solver is unconverged.  Solver can continue.
    UNCONVERGED,
    //! The status has not been checked.  Unknown state.
    UNCHECKED
  };

  enum CheckType {
    //! Check all status tests.
    COMPLETE,
    //! Only check critical tests and ignore optional or potentially costly tests.
    MINIMAL,
    //! Check is disabled.
    NONE
  };

  /** Pure virtual interface for a specifying a stopping criteria for the coupled system (strategy design pattern). */
  class StatusTest : public Teuchos::Describable,
		     public Teuchos::VerboseObject<pike::StatusTest> {

  public:

    /** \brief Test the stopping criteria.

	The test can (and should, if possible) be skipped if 
	checkType is "None".  If the test is skipped, then
	the status should be set to NOX::StatusTest::Unevaluated.
    */
    virtual pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE) = 0;
    
    //! Return the result of the most recent checkStatus() call.
    virtual pike::SolveStatus getStatus() const = 0;
    
    //! Reset the status test to reuse in another solve.
    virtual void reset() = 0;

  };
  
  /** \brief Prints the status test to the ostream
      \relates StatusTest
   */
  std::ostream& operator<<(std::ostream& os, const pike::StatusTest& test);

  std::string statusToString(const pike::SolveStatus& s);
  
}

#endif

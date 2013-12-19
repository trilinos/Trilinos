#ifndef PIKE_STATUS_TESTS_HPP
#define PIKE_STATUS_TESTS_HPP

#include <iostream>

namespace pike {

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

  class StatusTest {

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

    //! Output formatted description of stopping test to output stream.
    virtual std::ostream& print(std::ostream& stream) const = 0;

  };
  
  /** \brief Prints the status test to the ostream
      \relates StatusTest
   */
  std::ostream& operator<<(std::ostream& os, const pike::StatusTest& test);

}

#endif

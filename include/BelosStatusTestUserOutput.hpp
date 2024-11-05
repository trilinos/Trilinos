// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//

#ifndef BELOS_STATUS_TEST_USER_OUTPUT_HPP
#define BELOS_STATUS_TEST_USER_OUTPUT_HPP

/*!
  \file BelosStatusTestUserOutput.hpp
  \brief Special StatusTest for printing status tests in simple format for residuals.
*/

#include <vector>
#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

#include "BelosStatusTest.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestOutput.hpp"

namespace Belos {

  /*!
    \class StatusTestUserOutput
    \brief A special StatusTest for printing other status tests in a simple format.

  */
template <class ScalarType, class MV, class OP>
class StatusTestUserOutput : public StatusTestOutput<ScalarType,MV,OP> {

  typedef MultiVecTraits<ScalarType,MV> MVT;
  typedef Belos::StatusTestCombo<ScalarType,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestResNorm<ScalarType,MV,OP>  StatusTestResNorm_t;
  typedef Belos::StatusTestMaxIters<ScalarType,MV,OP>  StatusTestMaxIters_t;

 public:
  //! @name Constructors/destructors
  //@{

  /*! \brief Constructor
   *
   * The StatusTestUserOutput requires an OutputManager for printing the underlying StatusTest on
   * calls to checkStatus(), as well as an underlying StatusTest.
   *
   * The last two parameters, described below, in addition to the verbosity level of the OutputManager, control when printing is
   * called. When both the \c mod criterion and the \c printStates criterion are satisfied, the status test will be printed to the
   * OutputManager with ::MsgType of ::StatusTestDetails.
   *
   * @param[in] test A reference-counted pointer to an object that combines Belos::StatusTestMaxIters and a set of Belos::StatusTestResNorm
   * @param[in] mod A positive number describes how often the output should be printed. On every call to checkStatus(), an internal counter
   *                is incremented. Printing may only occur when this counter is congruent to zero modulo \c mod. Default: 1 (attempt to print on every call to checkStatus())
   * @param[in] printStates A combination of ::StatusType values for which the output may be printed. Default: ::Passed (attempt to print whenever checkStatus() will return ::Passed)
   *
   */
  StatusTestUserOutput(const Teuchos::RCP<OutputManager<ScalarType> > &printer,
      Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test,
      Teuchos::RCP<std::map<std::string,Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > > taggedTests,
      int mod = 1,
      int printStates = Passed)
    : printer_(printer),
      taggedTests_(taggedTests),
      state_(Undefined),
      headerPrinted_(false),
      stateTest_(printStates),
      modTest_(mod),
      lastNumIters_(-1),
      comboType_(0),
      blockSize_(1),
      currNumRHS_(0),
      currLSNum_(0),
      numLSDgts_(1),
      numIterDgts_(1)
    {
      // Set the input test.
      setChild(test);
    }

  //! Destructor
  virtual ~StatusTestUserOutput() {};
  //@}

  //! @name Status methods
  //@{
  /*! Check and return status of underlying StatusTest.

    This method calls checkStatus() on the StatusTest object passed in the constructor. If appropriate, the
    method will follow this call with a call to print() on the underlying object, using the OutputManager passed via the constructor
    with verbosity level ::StatusTestDetails.

    The internal counter will be incremented during this call, but only after
    performing the tests to decide whether or not to print the underlying
    StatusTest. This way, the very first call to checkStatus() following
    initialization or reset() will enable the underlying StatusTest to be
    printed, regardless of the mod parameter, as the current number of calls
    will be zero.

    If the specified Teuchos::RCP for the child class is Teuchos::null, then calling checkStatus() will result in a StatusTestError std::exception being thrown.

    \return ::StatusType indicating whether the underlying test passed or failed.
  */
  StatusType checkStatus( Iteration<ScalarType,MV,OP>* solver )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(iterTest_ == Teuchos::null,StatusTestError,"StatusTestUserOutput::checkStatus():  iteration test pointer is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(resTestVec_.size() == 0,StatusTestError,"StatusTestUserOutput::checkStatus():  residual test pointer is null.");
    state_ = test_->checkStatus(solver);

    // Update some information for the header, if it has not printed or the linear system has changed.
    LinearProblem<ScalarType,MV,OP> currProb = solver->getProblem();
    //if (!headerPrinted_ || currLSNum_ != currProb.getLSNumber()) {
    if (currLSNum_ != currProb.getLSNumber()) {
      currLSNum_ = currProb.getLSNumber();
      blockSize_ = solver->getBlockSize();
      currIdx_ = currProb.getLSIndex();
      currNumRHS_ = currIdx_.size();
      numLSDgts_ = (int)std::floor((double)MVT::GetNumberVecs(*(currProb.getRHS())))+1;
      numIterDgts_ = (int)std::floor(std::log10((double)iterTest_->getMaxIters()))+1;
    }
    // Print out current iteration information if it hasn't already been printed, or the status has changed
    if (((iterTest_->getNumIters() % modTest_ == 0) && (iterTest_->getNumIters()!=lastNumIters_)) || (state_ == Passed)) {
      lastNumIters_ = iterTest_->getNumIters();
      if ( (state_ & stateTest_) == state_) {
        if ( printer_->isVerbosity(StatusTestDetails) ) {
          print( printer_->stream(StatusTestDetails) );
        }
        else if ( printer_->isVerbosity(Debug) ) {
          print( printer_->stream(Debug) );
        }
      }
    }

    return state_;
  }

  //! Return the result of the most recent checkStatus call, or undefined if it has not been run.
  StatusType getStatus() const {
    return state_;
  }
  //@}


  //! @name Accessor methods
  //@{

  /*! \brief Set the output manager.
   */
  void setOutputManager(const Teuchos::RCP<OutputManager<ScalarType> > &printer) { printer_ = printer; }

  /*! \brief Set how often the child test is printed.
   */
  void setOutputFrequency(int mod) { modTest_ = mod; }

  /*! \brief Set child test, which must be a combination of a Belos::StatusTestMaxIters AND a combination of Belos::StatusTestResNorms and user-specified tests
   *
   *  \note This also resets the test status to ::Undefined.
   */
  void setChild(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test) {

    // First check to see if this test is a combination test
    Teuchos::RCP<StatusTestCombo_t> comboTest = Teuchos::rcp_dynamic_cast<StatusTestCombo_t>(test);
    TEUCHOS_TEST_FOR_EXCEPTION(comboTest == Teuchos::null,StatusTestError,"StatusTestUserOutput::setChild: The parameter \"test\" must be a Belos::StatusTestCombo.");
    std::vector<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > tmpVec = comboTest->getStatusTests();

    // Get the number of tests.
    int numTests = tmpVec.size();

    // Find the maximum iteration and residual tests
    for (int i=0; i<numTests; ++i) {

      // Check if this is a maximum iteration test.
      Teuchos::RCP<StatusTestMaxIters_t> tmpItrTest = Teuchos::rcp_dynamic_cast<StatusTestMaxIters_t>(tmpVec[i]);
      if (tmpItrTest != Teuchos::null) {
        iterTest_ = tmpItrTest;
        continue;
      }

      // Check if this is a single residual test
      // this should only be the case if there are no user-specified tests when using the hard-coded solver-
      // specific tests only
      Teuchos::RCP<StatusTestResNorm_t> tmpResTest = Teuchos::rcp_dynamic_cast<StatusTestResNorm_t>(tmpVec[i]);
      // If the residual status test is a single test, put in the vector
      if (tmpResTest != Teuchos::null) {
        resTestVec_.resize( 1 );
        resTestVec_[0] = tmpResTest;
        resTestNamesVec_.resize( 1 );
        resTestNamesVec_[0] = "IMPLICIT RES";
        continue;
      }

      // Check if the residual test is a combination of several StatusTestResNorm objects.
      // This should be the standard: we have a combo of the solver-specific tests and user-specific tests
      // The user-specific tests are probably a combo again.
      Teuchos::RCP<StatusTestCombo_t> tmpComboTest = Teuchos::rcp_dynamic_cast<StatusTestCombo_t>(tmpVec[i]);
      TEUCHOS_TEST_FOR_EXCEPTION(tmpComboTest == Teuchos::null,StatusTestError,"StatusTestUserOutput():  test must be Belos::StatusTest[MaxIters|ResNorm|Combo].");
      tmpVec = tmpComboTest->getStatusTests();
      comboType_ = tmpComboTest->getComboType();

      // Add only status tests which are not in the user-specified list of tagged status tests
      // More specifically: we want to add the implicit residual test
      typename std::map<std::string,Teuchos::RCP<StatusTest<ScalarType,MV,OP> > >::iterator it;
      for (size_t j=0; j<tmpVec.size(); ++j) {
        tmpResTest = Teuchos::rcp_dynamic_cast<StatusTestResNorm_t>(tmpVec[j]);

        if(tmpResTest == Teuchos::null) continue;

        bool bFound = false;
        for(it = taggedTests_->begin(); it != taggedTests_->end(); ++it) {
          if(tmpVec[j] == it->second) { bFound = true; break; }
        }
        if(!bFound) {
          resTestVec_.push_back(tmpResTest);
          resTestNamesVec_.push_back("IMPLICIT RES");
        }
      }

      // add all tagged tests (the ordering is by the Tag names in alphabetical ordering)
      for(it = taggedTests_->begin(); it != taggedTests_->end(); ++it) {
        resTestVec_.push_back(it->second);
        resTestNamesVec_.push_back(it->first);
      }
    }

    // Keep the pointer to the new test and reset the state to Undefined.
    test_ = test;
    state_ = Undefined;

  }

  //! \brief Get child test.
  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > getChild() const {
    return test_;
  }

  /*! \brief Set a short solver description for output clarity.
   */
  void setSolverDesc(const std::string& solverDesc) { solverDesc_ = solverDesc; }

  /*! \brief Set a short preconditioner description for output clarity.
   */
  void setPrecondDesc(const std::string& precondDesc) { precondDesc_ = precondDesc; }
  //@}


  //! @name Reset methods
  //@{
  /*! \brief Informs the status test that it should reset its internal configuration to the uninitialized state.
   *
   *  This resets the cached state to an ::Undefined state and calls reset() on the underlying test. It also
   *  resets the counter for the number of calls to checkStatus().
   */
  void reset() {
    state_ = Undefined;
    test_->reset();
    lastNumIters_ = -1;
    headerPrinted_ = false;
  }

  //! Informs the outputting status test that it should reset the number of calls to zero.
  /*! \note This outputting status test relies on the number of iterations performed, not number of times
      the status test has been called, so this method has no effect.
  */
  void resetNumCalls() {}

  //@}

  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const {
    std::string ind(indent,' ');
    std::string starLine(55,'*');
    std::string starFront(5,'*');

    std::ios_base::fmtflags osFlags(os.flags());

    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Print header if this is the first call to this output status test.
    if (!headerPrinted_) {
      os << std::endl << ind << starLine << std::endl;
      os << ind << starFront << " Belos Iterative Solver: " << solverDesc_ << std::endl;
      if (precondDesc_ != "")
        os << ind << starFront << " Preconditioner: " << precondDesc_ << std::endl;
      os << ind << starFront << " Maximum Iterations: " << iterTest_->getMaxIters() << std::endl;
      os << ind << starFront << " Block Size: " << blockSize_ << std::endl;
      os << ind << starFront << " Status tests: " << std::endl;
      test_->print(os,indent + 3);
      os << ind << starLine << std::endl;
      os.setf(std::ios_base::right, std::ios_base::adjustfield);
      std::string indheader( 7 + numIterDgts_, ' ' );
      os << indheader;
      for (int i=0; i<currNumRHS_; ++i) {
        if ( i > 0 && currIdx_[i]!=-1 ) {
          // Put in space where 'Iter :' is in the previous lines
          os << ind << indheader;
        }
        os << "[" << std::setw(numLSDgts_) << currIdx_[i]+1 << "] : ";
        for (size_t j=0; j<resTestVec_.size(); ++j) {
          os << std::setw(15) << resTestNamesVec_[j];
        }
        os << std::endl;
      }
      headerPrinted_ = true;
    }

    // Print out residuals for each residual test.
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::string ind2( 7 + numIterDgts_, ' ' );
    os << ind << "Iter " << std::setw(numIterDgts_) << iterTest_->getNumIters() << ", ";
    for (int i=0; i<currNumRHS_; ++i) {
      if ( i > 0 && currIdx_[i]!=-1 ) {
        // Put in space where 'Iter :' is in the previous lines
        os << ind << ind2;
      }
      os << "[" << std::setw(numLSDgts_) << currIdx_[i]+1 << "] : ";
      for (size_t j=0; j<resTestVec_.size(); ++j) {
        if ( resTestVec_[j]->getStatus() != Undefined && currIdx_[i]!=-1 ) {
          // distinguish between ResNormTest derived and others
          Teuchos::RCP<StatusTestResNorm_t> tempResTest = Teuchos::rcp_dynamic_cast<StatusTestResNorm_t>(resTestVec_[j]);
          if(tempResTest != Teuchos::null)
            os << std::setw(15) << (*tempResTest->getTestValue())[currIdx_[i]];
          else {
            if(resTestVec_[j]->getStatus() == Belos::Passed)
              os << std::setw(15) << "Passed";
            else if(resTestVec_[j]->getStatus() == Belos::Failed)
              os << std::setw(15) << "Failed";
            else os << std::setw(15) << "Undefined";
          }
        } else {
          os << std::setw(15) << "---";
        }
      }
      os << std::endl;
    }
    // reset os format
    os.flags(osFlags);
  }

  //@}

  private:
    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;

    // Overall status test.
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test_;

    // tagged tests
    Teuchos::RCP<std::map<std::string,Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > > taggedTests_;

    // Iteration test (as passed in).
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > iterTest_;

    //! Vector of residual status tests
    std::vector<Teuchos::RCP<StatusTest<ScalarType,MV,OP> > > resTestVec_;

    //! Name tags for status tests
    std::vector<std::string> resTestNamesVec_;

    std::string solverDesc_;
    std::string precondDesc_;
    std::vector<int> currIdx_;
    StatusType state_;
    mutable bool headerPrinted_;
    int stateTest_, modTest_;
    int lastNumIters_, comboType_;
    int blockSize_;
    int currNumRHS_, currLSNum_;
    int numLSDgts_, numIterDgts_;
};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_USER_OUTPUT_HPP */

// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//

#ifndef BELOS_STATUS_TEST_GENERAL_OUTPUT_HPP
#define BELOS_STATUS_TEST_GENERAL_OUTPUT_HPP

/*!
  \file BelosStatusTestGeneralOutput.hpp
  \brief Special StatusTest for printing any kind of status test.
*/


#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

#include "BelosStatusTest.hpp"
#include "BelosStatusTestOutput.hpp"
#include "BelosOutputManager.hpp"

namespace Belos {

  /*! 
    \class StatusTestGeneralOutput
    \brief A special StatusTest for printing other status tests. 
    
    StatusTestGeneralOutput is a wrapper around any another StatusTest that calls 
    StatusTest::print() on the underlying object on calls to StatusTestGeneralOutput::checkStatus().
    The frequency and occasion of the printing can be dictated according to some parameters passed to 
    StatusTestGeneralOutput::StatusTestGeneralOutput().
  */
template <class ScalarType, class MV, class OP>
class StatusTestGeneralOutput : public StatusTestOutput<ScalarType,MV,OP> {

 public:
  //! @name Constructors/destructors
  //@{ 

  /*! \brief Constructor
   *
   * The StatusTestGeneralOutput requires an OutputManager for printing the underlying StatusTest on
   * calls to checkStatus(), as well as an underlying StatusTest.
   *
   * The last two parameters, described below, in addition to the verbosity level of the OutputManager, control when printing is 
   * called. When both the \c mod criterion and the \c printStates criterion are satisfied, the status test will be printed to the 
   * OutputManager with ::MsgType of ::StatusTestDetails.
   *
   * @param[in] mod A positive number describes how often the output should be printed. On every call to checkStatus(), an internal counter
   *                is incremented. Printing may only occur when this counter is congruent to zero modulo \c mod. Default: 1 (attempt to print on every call to checkStatus())
   * @param[in] printStates A combination of ::StatusType values for which the output may be printed. Default: ::Passed (attempt to print whenever checkStatus() will return ::Passed)
   *
   */
  StatusTestGeneralOutput(const Teuchos::RCP<OutputManager<ScalarType> > &printer, 
			  Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test,
			  int mod = 1,
			  int printStates = Passed)
    : printer_(printer), 
      test_(test), 
      state_(Undefined), 
      stateTest_(printStates), 
      modTest_(mod), 
      numCalls_(0) 
  {}

  //! Destructor
  virtual ~StatusTestGeneralOutput() {};
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
  StatusType checkStatus( Iteration<ScalarType,MV,OP>* solver ) {
    TEUCHOS_TEST_FOR_EXCEPTION(test_ == Teuchos::null,StatusTestError,"StatusTestGeneralOutput::checkStatus(): child pointer is null.");
    state_ = test_->checkStatus(solver);

    if (numCalls_++ % modTest_ == 0) {
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

  /*! \brief Set child test.
   *
   *  \note This also resets the test status to ::Undefined.
   */
  void setChild(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test) {
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
    numCalls_ = 0;
  }

  //! Informs the outputting status test that it should reset the number of calls to zero.
  void resetNumCalls() { numCalls_ = 0; }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  void print(std::ostream& os, int indent = 0) const {
    std::string ind(indent,' ');
    os << std::endl << ind << "Belos::StatusTestGeneralOutput: ";
    switch (state_) {
    case Passed:
      os << "Passed" << std::endl;
      break;
    case Failed:
      os << "Failed" << std::endl;
      break;
    case Undefined:
      os << "Undefined" << std::endl;
      break;
    }
    os << ind << "  (Num calls,Mod test,State test): " << "(" << numCalls_ << ", " << modTest_ << ",";
    if (stateTest_ == 0) {
      os << " none)" << std::endl;
    }
    else {
      if ( stateTest_ & Passed ) os << " Passed";
      if ( stateTest_ & Failed ) os << " Failed";
      if ( stateTest_ & Undefined ) os << " Undefined";
      os << ")" << std::endl;
    }
    // print child, with extra indention
    test_->print(os,indent+3);
  }
 
  //@}

  private:
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test_;
    std::string solverDesc_;
    std::string precondDesc_;
    StatusType state_;
    int stateTest_;
    int modTest_;
    int numCalls_;
};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_OUTPUT_HPP */

//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//

#ifndef BELOS_STATUS_TEST_OUTPUT_HPP
#define BELOS_STATUS_TEST_OUTPUT_HPP

/*!
  \file BelosStatusTestOutput.hpp
  \brief Virtual base class for StatusTest that printing status tests.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

#include "BelosStatusTest.hpp"
#include "BelosOutputManager.hpp"


namespace Belos {

  /*! 
    \class StatusTestOutput
    \brief A virtual base class for StatusTest that print other status tests.     

    StatusTestOutput provides an interface for status tests that wrap another StatusTest.  These
    printing status tests can be generic, calling StatusTest::print() on the underlying object, or 
    can specifically require input status tests to be of a certain type. The frequency and occasion 
    of the printing can be dictated according to some parameters passed to 
    StatusTestOutput::StatusTestOutput().
  */
template <class ScalarType, class MV, class OP>
class StatusTestOutput : public virtual StatusTest<ScalarType,MV,OP> {

 public:
  //! @name Constructors/destructors
  //@{ 

  //! \brief Default constructor
  StatusTestOutput() {}

  /*! \brief Constructor
   *
   * The StatusTestOutput requires an OutputManager for printing the underlying StatusTest on
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
  StatusTestOutput(const Teuchos::RCP<OutputManager<ScalarType> > &printer, 
                   Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test,
                   int mod = 1,
                   int printStates = Passed)
    {}   

  //! Destructor
  virtual ~StatusTestOutput() {}
  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set the output manager.
   */ 
  virtual void setOutputManager(const Teuchos::RCP<OutputManager<ScalarType> > &printer) = 0;

  /*! \brief Set how often the child test is printed.
   */
  virtual void setOutputFrequency(int mod) = 0;

  /*! \brief Set child test.
   *
   *  \note This also resets the test status to ::Undefined.
   */
  virtual void setChild(Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test) = 0;

  //! \brief Get child test.
  virtual Teuchos::RCP<StatusTest<ScalarType,MV,OP> > getChild() const = 0;

  /*! \brief Set a short solver description for output clarity.
   */
  virtual void setSolverDesc(const std::string& solverDesc) = 0;

  /*! \brief Set a short preconditioner description for output clarity.
   */
  virtual void setPrecondDesc(const std::string& precondDesc) = 0;

  //@}


  //! @name Reset methods
  //@{ 

  //! Informs the outputting status test that it should reset the number of calls to zero.
  virtual void resetNumCalls() = 0;

  //@}

};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_OUTPUT_HPP */

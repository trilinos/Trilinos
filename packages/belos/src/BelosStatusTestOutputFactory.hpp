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

#ifndef BELOS_STATUS_TEST_OUTPUT_FACTORY_HPP
#define BELOS_STATUS_TEST_OUTPUT_FACTORY_HPP

/*!
  \file BelosStatusTestOutputFactory.hpp
  \brief A factory class for generating StatusTestOutput objects.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosStatusTestResNormOutput.hpp"
#include "BelosStatusTestGeneralOutput.hpp"


namespace Belos {

  /*! 
    \class StatusTestOutputFactory
    \brief A factory class for generating StatusTestOutput objects.

    StatusTestOutputFactory provides a generic interface for generating StatusTestOutput objects
    that any solver manager can use.  This factory removes the logic for selecting which StatusTestOutput
    class needs to be used from the solver managers. It also hides the generation of new StatusTestOutput
    classes from the solver managers. 
  */
template <class ScalarType, class MV, class OP>
class StatusTestOutputFactory {

 public:
  //! @name Constructors/destructors
  //@{ 

  /*! \brief Constructor
   *
   * The StatusTestOutputFactory generates a StatusTestOutput object that provides a particular style of
   * output, decided by the Belos::OutputType enumeration.
   *                
   * @param[in] outputStyle A ::OutputType value which defines the style of output requested by the user.
   */
  StatusTestOutputFactory( int outputStyle )
    : outputStyle_(outputStyle)
    {}

  //! Destructor
  virtual ~StatusTestOutputFactory() {}
  //@}


  //! @name Creation Methods
  //@{ 

  /*! \brief Create the StatusTestOutput object specified by the outputStyle
   *
   * The StatusTestOutput object requires an OutputManager for printing the underlying StatusTest on
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
   Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > create(const Teuchos::RCP<OutputManager<ScalarType> > &printer, 
                                                            Teuchos::RCP<StatusTest<ScalarType,MV,OP> > test,
                                                            int mod,
                                                            int printStates)
    {
      Teuchos::RCP<StatusTestOutput<ScalarType,MV,OP> > outputTest;

      switch( outputStyle_ ) {

      case General:
        if (mod > 0) {
          outputTest = Teuchos::rcp( new StatusTestGeneralOutput<ScalarType,MV,OP>( printer, test, mod, printStates ) );
        }
        else {
          outputTest = Teuchos::rcp( new StatusTestGeneralOutput<ScalarType,MV,OP>( printer, test, 1 ) );
        }
        break;
      case Brief:
        if (mod > 0) {
          outputTest = Teuchos::rcp( new StatusTestResNormOutput<ScalarType,MV,OP>( printer, test, mod, printStates ) );
        }
        else {
          outputTest = Teuchos::rcp( new StatusTestResNormOutput<ScalarType,MV,OP>( printer, test, 1 ) );
        }
        break;
      }

      return outputTest;
    }   

  //@}

 private:

  // Which type of StatusTestOutput class
  int outputStyle_;

  // Hide the default constructor and copy constructor
  StatusTestOutputFactory( void ) {}
  StatusTestOutputFactory( const StatusTestOutputFactory<ScalarType,MV,OP>& ) {}

};

} // end of Belos namespace

#endif /* BELOS_STATUS_TEST_OUTPUT_FACTORY_HPP */

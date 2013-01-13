// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_STATUSTEST_WRAPPER_H
#define LOCA_STATUSTEST_WRAPPER_H

#include "Teuchos_RCP.hpp"

#include "NOX_StatusTest_Generic.H" // base class

namespace LOCA {

  //! %LOCA %StatusTest classes
  namespace StatusTest {

    /*!
     * \brief A %wrapper class for wrapping a %NOX status test.
     */
    /*!
     * The %LOCA::StatusTest::Wrapper class provides a wrapper for %NOX
     * status tests to change the solver passed to the wrapped status test.
     * The solver passed through the checkStatus() method is wrapped via 
     * the LOCA::Solver::Wrapper class and then forwarded to the checkStatus()
     * method of the wrapped status test.  The purpose of this is to allow
     * status tests that use concrete group data to function correctly when
     * the group is stored in an extended continuation or bifurcation group.
     * (See LOCA::Solver::Wrapper for more details or the %LOCA status tests
     * page for examples on how to effectively use this class.)
     */
    class Wrapper : public NOX::StatusTest::Generic {

    public:

      //! Constructor
      Wrapper(const Teuchos::RCP<NOX::StatusTest::Generic>& s);

      //! Destructor
      virtual ~Wrapper();

      //! Calls checkStatus of underlying status test
      virtual NOX::StatusTest::StatusType
      checkStatus(const NOX::Solver::Generic& problem, 
		  NOX::StatusTest::CheckType checkType);

      //! Calls getStatus of underlying status test
      virtual NOX::StatusTest::StatusType getStatus() const;
      
      //! Calls print of underlying status test
      virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

      //! Returns underlying status test
      Teuchos::RCP<NOX::StatusTest::Generic>
      getUnderlyingStatusTest();

      //! Returns underlying status test
      Teuchos::RCP<const NOX::StatusTest::Generic>
      getUnderlyingStatusTest() const;

    private:

      //! Private to prohibit copying
      Wrapper(const Wrapper&);

      //! Private to prohibit copying
      Wrapper& operator=(const Wrapper&);
  
    protected:

      //! Stores pointer to underlying status test
      Teuchos::RCP<NOX::StatusTest::Generic> statusTestPtr;
      
    }; // class Wrapper

  } // namespace StatusTest
} // namespace LOCA


#endif

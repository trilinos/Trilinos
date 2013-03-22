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

#ifndef LOCA_STATUSTEST_COMBO_H
#define LOCA_STATUSTEST_COMBO_H


#include "LOCA_StatusTest_Abstract.H" // base class

#include "Teuchos_RCP.hpp"
#include "LOCA_GlobalData.H"


namespace LOCA {

namespace StatusTest {

/*!
  \brief Arbitrary combination of status tests.

  In the \c AND (see LOCA::StatusTest::Combo::ComboType) combination, the
  result is \c Unconverged (see LOCA::StatusTest::StatusType) if \e any of
  the tests is \c Unconverged. Otherwise, the result is equal to the
  result of the \e first test in the list that is either \c Converged
  or \c Failed. It is not recommended to mix \c Converged and \c
  Failed tests in an \c AND combination.

  In the \c OR combination, the result is \c Unconverged if \e all of
  the tests are \c Unconverged. Otherwise, it is the result of the \e
  first test in the list that is either \c Converged or \c
  Failed. Therefore, it will generally make sense to put the \c Failed
  -type tests at the end of the \c OR list.

  \note We call checkStatus on \e every convergence test, though some
  may be called with the LOCA::StatusTest::None option.
*/
class Combo : public Abstract {

public:

  /*! 
    \brief The test can be either the AND of all the component tests,
    or the OR of all the component tests.
  */
  enum ComboType {
    //! Logically "AND" together the results of the tests in this combination
    AND, 
    //! Logically "OR" together the results of the tests in this combination
    OR
  };

  //! Constructor.  Optional argument is the error stream for output.
  Combo( ComboType t,
         const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr = Teuchos::null );

  //! Constructor with a single test.
  Combo(ComboType t,
        const Teuchos::RCP<Abstract>& a,
        const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr = Teuchos::null );

  //! Constructor with two tests.
  Combo( ComboType t,
         const Teuchos::RCP<Abstract>& a,
         const Teuchos::RCP<Abstract>& b,
         const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr = Teuchos::null );

  //! Add another test to this combination. 
  /*!
    Calls isSafe() to determine if it is safe to add \c a to the combination.
  */
  virtual Combo& addStatusTest(const Teuchos::RCP<Abstract>& a);

  //! Destructor
  virtual ~Combo();

  /*!
    \brief Tests stopping criterion.
    
    See addOp() and orOp() for details.
  */
  virtual LOCA::StatusTest::StatusType
//  checkStatus( const LOCA::Stepper         & stepper,
  checkStatus( const LOCA::Abstract::Iterator         & stepper,
               LOCA::StatusTest::CheckType   checkType
             );

  virtual LOCA::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print( std::ostream& stream,
                          int indent = 0
                        ) const;
  
protected:

  //! Use this for checkStatus when this is an OR type combo. Updates LOCA::StatusTest::Combo::status.
  /*!  
    If there is a combination of LOCA::StatusTest::Failed and
    LOCA::StatusTest::Converged in the tests that are OR'd together,
    the value of status for this test is set to the status of the
    first test it encounters which is not LOCA::Status::Unconverged. The
    tests are evaluated in the order that they were added to the
    combination.

    \note We compute the status of all tests in the combination for
    the sake of completeness, even if we could determine the status of
    this combination test without that check.
  */
  virtual void
//  orOp( const LOCA::Stepper               & stepper,
  orOp( const LOCA::Abstract::Iterator               & stepper,
              LOCA::StatusTest::CheckType   checkType);

  //! Use this for checkStatus when this is an AND type combo. Updates LOCA::StatusTest::Combo::status.
  /*!  

    If any tests are LOCA::StatusTest::Unconverged, then the status of
    this test is LOCA::StatusTest::Unconverged.  If there is a
    combination of LOCA::StatusTest::Failed and
    LOCA::StatusTest::Converged in the tests that are AND'd together,
    the value of status for this test is set to the status of the
    first test it encounters.  The tests are evaluated in the
    order that they were added to the combination.

    \note We compute the status of all tests in the combination for
    the sake of completeness, even if we could determine the status of
    this combination test without that check.
  */
  virtual void
//  andOp( const LOCA::Stepper               & stepper,
  andOp( const LOCA::Abstract::Iterator               & stepper,
               LOCA::StatusTest::CheckType   checkType );

  /*! \brief Check whether or not it is safe to add \c a to this list
   of tests.

   This is necessary to avoid any infinite recursions 
   (i.e., a test cannot own a copy of itself).
  */
  bool isSafe(Abstract& a);

private:

  //! Type of test
  const ComboType type;

  //! Vector of generic status tests
  std::vector<Teuchos::RCP<Abstract> > tests;

  //! Status
  LOCA::StatusTest::StatusType status;

  //! Ostream used to print errors
  Teuchos::RCP<const LOCA::GlobalData> globalDataPtr_;

}; // class Combo

} // namespace Status
} // namespace LOCA


#endif

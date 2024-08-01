// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
